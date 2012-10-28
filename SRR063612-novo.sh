#!/usr/bin
# Mapping with Novocraft

echo ===============[`date`]==============

# Align with novoalign
novoalign -d ~/Data/Data/Rice/Os-Nipponbare-Reference-IRGSP-1.0/novo-IRGSP-ref -f Cultivar/Japonica/ARO/IRGC9060/SRR063612_1.fastq Cultivar/Japonica/ARO/IRGC9060/SRR063612_2.fastq -i PE 468,33 -F SLXFQ -r None -o SoftClip -o SAM > Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.sam

# Convert SAM to BAM
samtools view -bS Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.sam -o Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.bam

# Sort
samtools sort Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.bam Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.sorted 

# Call variant candidates
samtools mpileup -ugf ~/Data/Data/Rice/Os-Nipponbare-Reference-IRGSP-1.0/IRGSP-1.0_genome.fasta Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.sorted.bam | bcftools view -bvcg - > Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.var.raw.bcf

bcftools view Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.var.raw.bcf | vcfutils.pl varFilter -D 100 > Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.var.flt.vcf 

echo `date` ':all processes completed!'
echo

# do some cleaning
rm Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.sam
rm Cultivar/Japonica/ARO/IRGC9060/SRR063612.novo.pe.bam

echo ---------------------------------------

