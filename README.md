# Computational Human Genomics
This project is part of the course in Computational Human Genomics, held by Dr.Yari Ciani and Prof. Francesca Demichelis, during the A.A. 2022-2023. 
This project aims to apply an example of human genomic workflow to a patient, characterizing both tumor and control samples. Thanks to the applied workflow we were able to characterize germline and somatic variants, determine the ancestry of the patient, and study the tumor ploidy and purity.

Project developed by: 
  * Andrea Tonina  [@iamandreatonina](https://github.com/iamandreatonina)
  * Gloria Lugoboni [@GloriaLugoboni](https://github.com/GloriaLugoboni)
  * Lorenzo Santarelli [@Lor-Santa](https://github.com/Lor-Santa)

[Report]()

---

In particular, the project is aimed at answering 10 specific tasks:
 * Explore statistics about the raw aligned reads contained in the two BAM files;
 * Perform realignment and recalibration on those;
 * Identify and annotate heterozygous SNPs;
 * Determine the ancestry of the patients;
 * Identify somatic copy number variants;
 * Identify somatic point mutations;
 * Determine how DNA Repair Genes have been impacted by germline CNVs and SNPs;
 * Determine which DNA repair genes overlap both germline heterozygous copy-number deletions and somatic point mutations;
 * Determine tumor purity and ploidy;
 * Determine the similarity of Tumor and Control samples.

---

Tools utilized: 
- [X] Samtools
- [X] GATK 
- [X] Picard
- [X] BCFTOOLS
- [X] SnpEff
- [X] Varscan
- [X] EthSEQ
- [X] CBS
- [X] DNACopy
- [X] CLONET
- [X] TPES
- [X] SPIA
