
#atsisiunciame duomenis

#SRR11647648
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/048/SRR11647648/SRR11647648_1.fastq.gz &
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/048/SRR11647648/SRR11647648_2.fastq.gz &

#SRR11647649
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/049/SRR11647649/SRR11647649_1.fastq.gz &
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/049/SRR11647649/SRR11647649_2.fastq.gz &

#SRR11647658
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/058/SRR11647658/SRR11647658_1.fastq.gz &
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/058/SRR11647658/SRR11647658_2.fastq.gz &

#SRR11647659
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/059/SRR11647659/SRR11647659_1.fastq.gz &
wget -P raw_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/059/SRR11647659/SRR11647659_2.fastq.gz

ls -lh raw_data/*.fastq.gz
gzip -t raw_data/*.fastq.gz && echo "all ok"


#QC
mkdir fastqc_res
fastqc -t 8 raw_data/*.fastq.gz -o fastqc_res

mkidir multiqc_res
multiqc fastqc_res -o multiqc_res

#trimming:
mkdir ~/homework2/trimmed_data
cd ~/homework2/raw_data
for sample in SRR11647648 SRR11647649 SRR11647658 SRR11647659; do
    trim_galore --paired --cores 8 -o ~/homework2/trimmed_data/ ${sample}_1.fastq.gz ${sample}_2.fastq.gz
done
cat trimmed_data/*_trimming_report.txt | grep "Reads with adapters"


mkdir trimmed_fastqc
fastqc -t 8 trimmed_data/*.fq.gz -o trimmed_fastqc
#trimming stats
grep -E "reads processed|Reads with adapters|Quality-trimmed|Total written" ~/homework2/trimmed_data/*_trimming_report.txt > ~/homework2/trimmed_data/trimming_summary.txt

#kazkodel su viena seka nepadare
fastqc -t 8 trimmed_data/SRR11647648_1_val_1.fq.gz -o trimmed_fastqc

multiqc trimmed_fastqc -o multiqc_res


#mapping
mkdir mapped
cd trimmed_data

for sample in SRR11647648 SRR11647649 SRR11647658 SRR11647659; do
    bsmap -a ${sample}_1_val_1.fq.gz -b ${sample}_2_val_2.fq.gz \
    -d ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    -o ~/homework2/mapped/${sample}.bam \
    -p 6
done

ls -lh ~/homework2/mapped/

# kazkodel vieno nesumapino iki galo.....
rm ~/homework2/mapped/SRR11647659.bam
tmux new -s mapping
conda activate bsmap_env
cd ~/homework2/trimmed_data
bsmap -a SRR11647659_1_val_1.fq.gz -b SRR11647659_2_val_2.fq.gz \
    -d ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    -o ~/homework2/mapped/SRR11647659.bam \
    -p 6


#marking duplicates
for SAMPLE in SRR11647648 SRR11647649 SRR11647658 SRR11647659; do
    samtools collate -@ 6 -O -u mapped/${SAMPLE}.bam | \
    samtools fixmate -@ 6 -m -u - - | \
    samtools sort -@ 6 -u - | \
    samtools markdup -@ 6 - mapped/${SAMPLE}_markdup.bam
    echo "${SAMPLE} duplicates:" >> mapped/qc_stats.txt
    samtools view -c -f 1024 mapped/${SAMPLE}_markdup.bam >> mapped/qc_stats.txt
    echo "${SAMPLE} flagstat:" >> mapped/qc_stats.txt
    samtools flagstat mapped/${SAMPLE}_markdup.bam >> mapped/qc_stats.txt
    echo "---" >> mapped/qc_stats.txt
done

grep "Total reads processed" ~/homework2/trimmed_data/*_trimming_report.txt
grep "properly paired" ~/homework2/mapped/qc_stats.txt

#indexing mapped data
cd ~/homework2/mapped
for sample in SRR11647648 SRR11647649 SRR11647658 SRR11647659; do
    samtools sort -@ 6 ${sample}_markdup.bam -o ${sample}_markdup_sorted.bam
    samtools index ${sample}_markdup_sorted.bam
done

ls -lh ~/homework2/mapped/*_markdup_sorted.bam
ls ~/homework2/mapped/*_markdup_sorted.bam.bai

conda create -n methyldackel_env
conda activate methyldackel_env
conda install -c bioconda methyldackel



#mbias
mkdir ~/homework2/mbias_20chr
conda activate methyldackel_env

samtools view -H ~/homework2/mapped/SRR11647648_markdup_sorted.bam | grep "^@SQ"

for sample in SRR11647648 SRR11647649 SRR11647658 SRR11647659; do
    MethylDackel mbias \
        -r 20 \
        ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        ~/homework2/mapped/${sample}_markdup_sorted.bam \
        ~/homework2/mbias_20chr/${sample}_mbias_20chr
done


#methylation calling
mkdir ~/homework2/methylation_chr20

MethylDackel extract --OT 0,0,2,120 --OB 2,0,6,124 -r 20 \
    ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    ~/homework2/mapped/SRR11647648_markdup_sorted.bam \
    -o ~/homework2/methylation_chr20/SRR11647648_methylation

MethylDackel extract --OT 2,0,2,112 --OB 2,124,15,124 -r 20 \
    ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    ~/homework2/mapped/SRR11647649_markdup_sorted.bam \
    -o ~/homework2/methylation_chr20/SRR11647649_methylation

MethylDackel extract --OT 3,0,2,114 --OB 2,123,6,124 -r 20 \
    ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    ~/homework2/mapped/SRR11647658_markdup_sorted.bam \
    -o ~/homework2/methylation_chr20/SRR11647658_methylation

MethylDackel extract --OT 2,0,2,120 --OB 2,123,6,124 -r 20 \
    ~/homework2/references/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    ~/homework2/mapped/SRR11647659_markdup_sorted.bam \
    -o ~/homework2/methylation_chr20/SRR11647659_methylation

ls -lh ~/homework2/methylation_chr20/

head -5 ~/homework2/methylation_chr20/SRR11647648_methylation_CpG.bedGraph