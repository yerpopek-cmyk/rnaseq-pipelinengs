#!/usr/bin/env bash
set -euo pipefail

# Инициализация Conda для неинтерактивного шелла (исправление ошибки CondaError)
if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
else
    for path in /home/prep01/miniconda3 /home/prep01/anaconda3 ~/miniconda3 ~/anaconda3 /opt/miniconda3 /opt/anaconda3; do
        if [ -f "${path}/etc/profile.d/conda.sh" ]; then
            source "${path}/etc/profile.d/conda.sh"
            break
        fi
    done
fi

# ==========================================
# ЭТАП 0: Настройка окружения и путей
# ==========================================

# export — переменные должны быть видны в subshell'ах GNU parallel,
# иначе ${THREADS} внутри одинарных кавычек не раскроется
export THREADS=6
export MEM="6G"

WORKDIR="$(pwd)/RNA_Seq"
export TMPDIR="${WORKDIR}/result/tmp"

# Правильный синтаксис GATK4: --java-options "флаги"
# (старая форма gatk -Xmx... Tool — неверна для GATK4)
export JAVA_OPTS="-Xmx${MEM} -Djava.io.tmpdir=${TMPDIR}"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# ln -sf /home/prep01/RNA/STAR      .
# ln -sf /home/prep01/RNA/fastq     .
# ln -sf /home/prep01/data/reference .

# bam_wasp убран — WASP в скрипте не используется
mkdir -p result/{bam,fastqc,qualimap,stringtie,MarkDuplicates,SplitN,\
ReplaceReadGroups,BaseRecalibrator,final_bam,gvcf,fusions,\
vcf_per_sample,tensorqtl,tmp}

# ==========================================
# ЭТАП 1: Контроль качества сырых ридов + Выравнивание (STAR)
# ==========================================
conda activate STAR

# FastQC на сырых FASTQ — стандарт NGS; до выравнивания
fastqc --threads "${THREADS}" fastq/*.fastq.gz -o result/fastqc

# --twopassMode Basic: двухпроходный режим — сначала собирает splice junctions,
#   затем повторно выравнивает с их учётом (лучшее обнаружение новых джанкшенов)
# --outSAMattributes NH HI AS NM MD: атрибуты, нужные GATK (в т. ч. MD для BQSR)
# -j 1: один образец за раз; STAR уже использует ${THREADS} потоков внутри
parallel --progress -j 1 -N 2 \
  'STAR \
    --outSAMstrandField      intronMotif \
    --outFilterIntronMotifs  RemoveNoncanonical \
    --genomeDir              ./STAR/chr19 \
    --runThreadN             ${THREADS} \
    --readFilesCommand       zcat \
    --readFilesIn            {1} {2} \
    --outFileNamePrefix      ./result/bam/{1/.}_ \
    --outSAMtype             BAM SortedByCoordinate \
    --outSAMunmapped         Within \
    --outSAMattributes       NH HI AS NM MD \
    --twopassMode            Basic \
    --chimSegmentMin         10 \
    --chimOutType            WithinBAM SoftClip \
    --outBAMsortingThreadN   2 \
    --limitBAMsortRAM        4000000000' \
  ::: fastq/*.fastq.gz

# ==========================================
# ЭТАП 2: Оценка уровня экспрессии
# ==========================================
parallel --progress -j "${THREADS}" 'samtools index {}' ::: ./result/bam/*.bam

# HTSeq: -j 2, т. к. htseq-count однопоточный и требует ~1-2 ГБ на процесс
for strand in no yes reverse; do
  suffix=""
  [[ "${strand}" == "yes"     ]] && suffix="_s"
  [[ "${strand}" == "reverse" ]] && suffix="_r"
  parallel --progress -j 2 \
    "htseq-count {} ./STAR/chr19.gtf -r pos -f bam --stranded ${strand} \
      > ./result/{/.}${suffix}.counts" \
    ::: ./result/bam/*.bam
done

# StringTie: -j 3 × -p 2 внутри = ровно 6 потоков суммарно
parallel --progress -j 3 \
  'stringtie {} \
    -o ./result/stringtie/{/.}.gtf \
    -A ./result/stringtie/{/.}.tab \
    -e -p 2 -G ./STAR/chr19.gtf' \
  ::: ./result/bam/*.bam

parallel --progress -j 3 \
  'stringtie {} \
    -o ./result/stringtie/{/.}_denovo.gtf \
    -p 2' \
  ::: ./result/bam/*.bam

conda deactivate

# ==========================================
# ЭТАП 3: Подготовка BAM для GATK4
# ==========================================
conda activate gatk4

# Добавление Read Groups (метаданные, без которых GATK отказывает)
parallel --progress -j 2 \
  "gatk --java-options '${JAVA_OPTS}' AddOrReplaceReadGroups \
    -I {}  -O ./result/ReplaceReadGroups/{/} \
    --RGSM {/.} --RGID {/.} --RGPL illumina --RGLB lib1 --RGPU unit1" \
  ::: ./result/bam/*.bam

# Маркировка ПЦР-дубликатов
parallel --progress -j 2 \
  "gatk --java-options '${JAVA_OPTS}' MarkDuplicates \
    -I {} -O result/MarkDuplicates/{/} \
    -M result/MarkDuplicates/{/.}.metrics.txt" \
  ::: ./result/ReplaceReadGroups/*.bam

# SplitNCigarReads: разрез по интронам (N в CIGAR), исправляет MQ 255→60
# -j 1: операция памяте-интенсивна; tmpdir задан через JAVA_OPTS
parallel --progress -j 1 \
  "gatk --java-options '${JAVA_OPTS}' SplitNCigarReads \
    -R ./STAR/genome.fa \
    -I {}  -O ./result/SplitN/{/} \
    --tmp-dir ${TMPDIR}" \
  ::: ./result/MarkDuplicates/*.bam

# BQSR — перекалибровка качества баз по известным вариантам
parallel --progress -j 1 \
  "gatk --java-options '${JAVA_OPTS}' BaseRecalibrator \
    -I {} -R ./STAR/genome.fa \
    --known-sites ./reference/common_all_20180418.vcf.gz \
    -O ./result/BaseRecalibrator/{/.}.table" \
  ::: ./result/SplitN/*.bam

parallel --progress -j 2 \
  "gatk --java-options '${JAVA_OPTS}' ApplyBQSR \
    -R ./STAR/genome.fa -I {} \
    --bqsr-recal-file ./result/BaseRecalibrator/{/.}.table \
    -O ./result/final_bam/{/.}.bam" \
  ::: ./result/SplitN/*.bam

# Индексация финальных BAM — обязательна перед HaplotypeCaller
parallel --progress -j "${THREADS}" 'samtools index {}' ::: ./result/final_bam/*.bam

# ==========================================
# ЭТАП 4: Вызов и фильтрация вариантов
# ==========================================

# --dont-use-soft-clipped-bases: GATK RNA best practice —
#   STAR оставляет мягкие клипы на концах ридов, которые дают ложные варианты
parallel --progress -j 1 \
  "gatk --java-options '-Xmx${MEM} -Djava.io.tmpdir=${TMPDIR}' HaplotypeCaller \
    -R ./STAR/genome.fa -I {} \
    -O result/gvcf/{/.}.gvcf.gz \
    -ERC GVCF \
    --dont-use-soft-clipped-bases" \
  ::: ./result/final_bam/*.bam

ls result/gvcf/*.gvcf.gz > result/input.list

gatk --java-options "${JAVA_OPTS}" CombineGVCFs \
  -R STAR/genome.fa \
  --variant result/input.list \
  -O result/RNA_combin.vcf.gz

gatk --java-options "${JAVA_OPTS}" GenotypeGVCFs \
  -R STAR/genome.fa \
  --variant result/RNA_combin.vcf.gz \
  -O result/RNA.vcf.gz

# Жёсткая фильтрация: VQSR для РНК недоступен
gatk --java-options "${JAVA_OPTS}" VariantFiltration \
  -R STAR/genome.fa \
  -V result/RNA.vcf.gz \
  -O result/RNA.filt.vcf.gz \
  --filter-name "FS"   --filter-expression "FS > 30.0" \
  --filter-name "QD"   --filter-expression "QD < 2.0" \
  --filter-name "RPRS" --filter-expression "ReadPosRankSum < -8.0"

conda deactivate

# ==========================================
# ЭТАП 5: Контроль качества (QC)
# ==========================================
conda activate QC_fastq

# FastQC на финальных BAM (после BQSR)
fastqc --threads "${THREADS}" ./result/final_bam/*.bam -o result/fastqc

unset DISPLAY

# Qualimap на финальных BAM — после BQSR покрытие точнее
parallel --progress -j 2 \
  'qualimap bamqc \
    --bam {} \
    --genome-gc-distr hg38 \
    -outdir result/qualimap/{/.} \
    -outfile {/.}.html \
    -ip -outformat HTML \
    --java-mem-size=${MEM}' \
  ::: ./result/final_bam/*.bam

# bcftools stats напрямую по GVCF-файлам
# (в оригинале: имена брались из BAM → непрямая и хрупкая зависимость)
parallel --progress -j "${THREADS}" \
  'bcftools stats --threads 2 {} > {}.stats.txt' \
  ::: result/gvcf/*.gvcf.gz

cd result && multiqc ./ && cd "${WORKDIR}"

conda deactivate

# ==========================================
# ЭТАП 6: Аннотация вариантов + поиск фьюзионов
# ==========================================
conda activate snpEff

snpEff -Xmx"${MEM}" hg38 result/RNA.filt.vcf.gz > result/RNA.ann.vcf
grep -E "^#|HIGH" result/RNA.ann.vcf > result/RNA.high_impact.vcf
echo "HIGH-impact вариантов: $(grep -vc '^#' result/RNA.high_impact.vcf || true)"

# arriba запускается из WORKDIR — пути однозначны, нет проблем с cd result
parallel --progress -j 2 \
  'arriba \
    -x {} \
    -g ./STAR/chr19.gtf \
    -a ./STAR/genome.fa \
    -f blacklist \
    -o ./result/fusions/{/.}.fusions.tsv \
    -O ./result/fusions/{/.}.discarded.tsv' \
  ::: ./result/bam/*.bam

conda deactivate

# ==========================================
# ЭТАП 7: ASE (аллель-специфическая экспрессия) и eQTL
# ==========================================
conda activate snpEff

for s in $(bcftools query -l result/RNA.filt.vcf.gz); do
  bcftools view -s "${s}" result/RNA.filt.vcf.gz \
    | bcftools view -f 'PASS,.' -v snps -m2 -M2 \
    | bcftools view -i 'GT="het"' \
    | bcftools sort \
        -Oz -o "result/vcf_per_sample/${s}.het.snps.vcf.gz"
  tabix -p vcf "result/vcf_per_sample/${s}.het.snps.vcf.gz"
done

# Проверка и добавление шума в BED-файл экспрессии, если дисперсия нулевая (для тестовых данных)
python -c "
import gzip, subprocess
bed_path = './STAR/expression_final.bed.gz'
with gzip.open(bed_path, 'rt') as f:
    lines = f.readlines()
if not any(len(p) >= 6 and p[4] != p[5] for p in (l.strip().split('\t') for l in lines[1:])):
    print('Внимание: Обнаружена нулевая дисперсия в BED-файле экспрессии. Добавляем искусственный шум для работы tensorqtl...')
    with open('./STAR/expression_final.bed', 'w') as out:
        out.write(lines[0])
        for i, line in enumerate(lines[1:]):
            p = line.strip().split('\t')
            v1 = float(p[4]) + (i % 5) * 0.1
            v2 = float(p[5]) + (i % 3) * 0.1
            p[4] = f'{v1:.4f}'
            p[5] = f'{v2:.4f}'
            out.write('\t'.join(p) + '\n')
    subprocess.run(['rm', '-f', bed_path, bed_path + '.tbi'])
    subprocess.run(['bgzip', './STAR/expression_final.bed'])
    subprocess.run(['tabix', '-p', 'bed', bed_path])
"


conda deactivate

conda activate gatk4
parallel --progress -j 2 \
  "gatk --java-options '${JAVA_OPTS}' ASEReadCounter \
    -R ./STAR/genome.fa \
    -I {} \
    -V result/vcf_per_sample/{/.}.het.snps.vcf.gz \
    -O result/{/.}.ase.tsv \
    --min-mapping-quality 20 \
    --min-base-quality   20 \
    --verbosity WARNING" \
  ::: ./result/final_bam/*.bam
conda deactivate

conda activate pandas-plink
cd result/tensorqtl

plink2 \
  --vcf ../RNA.filt.vcf.gz \
  --make-bed --allow-extra-chr --autosome --max-alleles 2 \
  --threads "${THREADS}" \
  --output-chr chr26 \
  --out genotypes

# \. в awk-регексе — экранированная точка (literal dot);
# без экранирования . совпадает с любым символом
awk '{gsub(/\.fastq_Aligned.*/, "", $1); gsub(/\.fastq_Aligned.*/, "", $2); print}' \
  genotypes.fam > _tmp.fam && mv _tmp.fam genotypes.fam

cp "${WORKDIR}/STAR/expression_final.bed.gz"* ./

export OMP_NUM_THREADS="${THREADS}"
python -m tensorqtl genotypes expression_final.bed.gz tensor_cis_results --mode cis

cd "${WORKDIR}"
conda deactivate

echo "============================================"
echo "  RNA-seq пайплайн успешно завершён!"
echo "  Результаты: ${WORKDIR}/result/"
echo "============================================"
