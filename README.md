This documentation accompanies the manuscript “Glacial refugia, post-glacial dynamics and hybrid zones of Pinaceae in Eurasia captured from sedimentary ancient DNA” by Stefano Meucci, Kathleen R. Stoof-Leichsenring, Yanrong Zhang, Andrei A. Andreev, Konstantin V. Krutovsky, Boris K. Biskaborn, Inger G. Alsos, Laura Parducci, Kevin Nota, Sandra Garcés-Pastor, Petr Kuneš, Walter Finsinger, Eleonora Cagliero, Daniel Vondrák, Jaroslav Kukla, Darrell Kaufman, Barbara Wohlfarth, Heikki Seppä, and Ulrike Herzschuh. The associated pipeline reflects the workflow presented in the manuscript and includes a series of Bash and R scripts for bioinformatic preprocessing and data analysis. The content and specific purpose of each script are outlined below.

Script: 01_FASTQC_FASTP.sl
Notes: Script for processing multiple samples of a single sediment core (e.g., Bolshoe Toko), to be adapted for each individual core.
Demultiplexed and adapter-trimmed FASTQ files, as provided by the sequencing company, were quality-checked using FASTQC (Andrews 2015) both before and after deduplication and trimming. Initial deduplication was performed with CLUMPIFY (v. 39.01) using the “clumpify.sh” script (Bushnell 2014), followed by trimming with FASTP (v. 0.23.2) (S. Chen et al. 2018) using the parameters: --merge, --length_required=30, --overlap_len_require=5, --correction, --low_complexity_filter, --cut_front, --cut_tail, --cut_window_size=4, and --cut_mean_quality=10. After merging, a second deduplication was carried out on both merged and unmerged reads using the “dedupe.sh” script from CLUMPIFY (v. 39.01) (Bushnell 2014).

Script: 02_KRAKEN_KRONA.sl
Notes: Script to process all sediment cores simultaneously.
Chloroplast-enriched DNA sequence reads were analyzed and classified with KRAKEN2 (v. 2.1.2) (Wood, Lu, and Langmead 2019) using a confidence threshold of 0.5 against the plastid database (RefSeq plastid release of NCBI (O’Leary et al. 2016)), downloaded in September 2023, with the addition of the newly assembled cp genome of P. obovata as described in the section 2.1. The database was built using KRAKEN2’s standard parameters (k-mer length of 35 bp and minimizer length of 31 bp). 

Script: 03_EXTRACT_PINACEAE_READS.sl
Notes: Script for processing multiple samples of a single sediment core (e.g., Bolshoe Toko), to be adapted for each individual core.
Merged and unmerged reads classified at the family level to Pinaceae were extracted using KRAKENTOOLS extract_kraken_reads.py (v. 1.2) (Lu et al. 2022).

Script: 04.1_READ_COUNT_PINACEAE.sl
Notes: Script to process all sediment cores simultaneously.
Script: 4.2_PLOT_READ_COUNT.R
Notes: The R script includes the plotting of both Pinaceae-assigned reads and genus-specific reads (therefore this same script is also provided below in the corresponding section).
The extracted Pinaceae reads were counted and plotted for each sample and site using the R (R Core Team 2022) package tidyverse (Wickham et al. 2019) (Appendix S1: Figures S1–S19; Appendix S2: Table S6).

Script: 05_MAPPING_TO_CONCATENATE_GENOME.sl
Notes: Script for processing multiple samples of a single sediment core (e.g., Bolshoe Toko), to be adapted for each individual core.
Reads classified as Pinaceae were aligned against a concatenated genome composed of the P. obovata (this study, under submission), P. pumila (NCBI GenBank: NC_041108.1), L. gmelinii (NCBI GenBank: MK468637.1) and A. sibirica (NCBI GenBank: NC_035067.1) cp reference genomes (Appendix S2: Table S3A) using BWA aln (v. 0.7.17) (H. Li and Durbin 2009)  with the parameters -l 1024 -o 2 -n 0.01, as recommended for ancient DNA read mapping by Oliva et al. (Oliva et al. 2021). Mapping quality and length filtering were performed using BAMTOOLS filter (Barnett et al. 2011), with the parameters:  -mapQuality '>=25' -length '>=35'. Further processing of the alignment files, including conversion, sorting and indexing was conducted using SAMTOOLS (v. 1.16.1) (H. Li et al. 2009). A third round of deduplication was performed, as deduplication is more effective on aligned reads. Deduplication for both merged and unmerged reads was performed using PICARD MarkDuplicates (v.3.1.0) (Broad Institute 2019). The alignments from the merged and unmerged reads were then combined into a single BAM file using SAMTOOLS merge (H. Li et al. 2009).

Script: 06_EXTRACT_GENUS_SPECIFIC_READS.sl
Notes: Script to process all sediment cores simultaneously.
Mapped reads to the concatenated genome were subsequently extracted into genus-specific BAM files using SAMTOOLS view (v. 1.16.1) (H. Li et al. 2009), by specifying the base pair region range (first and last base pair number) of the four cp reference genomes. 

Script: 07.1_CALCULATE_DEPTH_COVERAGE
Script: 07.2_REARRANGE_DEPTH_COVERAGE_RESULTS.R
Notes: Script to process all sediment cores simultaneously.
The breadth (bases covered by at least one read) and depth of coverage of bait regions were calculated using SAMTOOLS depth, followed by custom AWK scripts to compute the mean depth and breadth of coverage for each sample (Appendix S2: Table S6). 

Script: 08_READ_COUNT_GENUS_SPECIFIC_READS.sl
Script: 03.3_PLOT_READ_COUNT.R
Notes: Script to process all sediment cores simultaneously. The script includes the plotting of both Pinaceae-assigned reads and genus-specific reads (therefore this same script was already provided in the corresponding section 3).
Genus-specific reads were counted from the genus-specific BAM files and plotted for each sample and site using the R (R Core Team 2022) package tidyverse (Wickham et al. 2019) (Appendix S1: Figures S1–S19; Appendix S2: Table S6). 

Script: 09_GENUS_REFERENCE_MAPPING_SNP_CALL.sl
Notes: One script for all four genera, the references used as “REF” and “IN” are listed in Appendix S2: Table S3C and S3D, respectively.
To identify species-specific SNPs of each Pinaceae genus, available cp reference genomes from species expected in the defined geographical regions (Appendix S2: Table S3D) were mapped to cp reference genomes from species that are phylogenetically and geographically distant from those expected at our sites using BWA mem (v. 0.7.17) (H. Li and Durbin 2009) (e.g., Picea pungens, mostly located in central US, was selected due to its significant geographic distance and its different phylogenetic clade from the species likely present at our study sites; Lockwood et al., 2013) (Appendix S2: Table S3C). Variations were called using FREEBAYES (v. 1.3.6) (Garrison and Marth 2012) with parameter settings --min-alternate-count, 1 --min-alternate-fraction 0, --haplotype-length 0, and --pooled- continuous. 

Script: 12_SNP-BASED_SPECIES_IDENTIFICATION.R
Notes: This R script (Script 12) performs SNP-based species identification by comparing species-specific SNPs ascertained in Script 9 with the SNPs called from sediment samples in Script 11. Although it is executed as Script 12, it is included also here to align with the sequential structure of the Materials and Methods pipeline and is reiterated below as the twelfth step.
From the resulting SNP dataset, variant positions, reference alleles, and alternative alleles were extracted using R (R Core Team 2022). Ascertained species-specific SNPs were labeled accordingly (e.g., Picea obovata) (Appendix S2: Table S7). If multiple reference genomes from the same species were present, their SNPs were merged under the species name, ensuring that only species-level distinctions were maintained. When multiple reference species shared the same alternative allele at a given position, the SNP was excluded from the species determination analysis, as it was not phylogenetically informative. The number of species-specific SNPs decreased as more chloroplast reference genomes were included, due to the increasing number of common SNPs that were subsequently excluded. Therefore, the selection of reference genomes was restricted to species expected within the defined geographical regions (Appendix S2: Table S3D).

Script: 10_MAPPING_GENUS-SPECIFIC_READS.sl
Notes: Script for processing multiple samples of a single sediment core, to be adapted for each individual core.
The genus-specific reads from sediment samples were mapped to the same cp reference genomes from the phylogenetically and geographically distant species (Appendix S2: Table S3C) using BWA aln (v. 0.7.17) (H. Li and Durbin 2009) with the parameters -l 1024 -o 2 -n 0.01. Alignment files were processed with SAMTOOLS (v. 1.16.1) (H. Li and Durbin 2009) for conversion, sorting, and indexing. Read groups were added with PICARD AddOrReplaceReadGroups (v.3.1.0) (Broad Institute 2019). 

Script: 11.1_SNP_CALL_GENUS-SPECIFIC_READS.sl
Script: 11.2_RENAME_SAMPLES_IN_VCF_FILE.R
Notes: Script 11.1 performs SNP calling on genus-specific reads mapped to chloroplast reference genomes from phylogenetically and geographically distant species. This is done separately for each of the four genera, but includes all sediment cores simultaneously. Script 11.2 is used to rename samples in the resulting genus-specific VCF files.
Genus-specific SNPs were jointly called for all sediment samples using FREEBAYES (v. 1.3.6) (Garrison and Marth 2012) with the options --pooled-continuous and --min-base-quality 10. The resulting VCF files were processed with PLINK (v. 1.90) (Chang et al. 2015) to remove samples with >90% missing data, while VCFTOOLS (v. 0.1.16) (Danecek et al. 2011) was used to filter out indels, sites with >70% missing data, and minor allele frequencies below 1%. 

Script: 12_SNP-BASED_SPECIES_IDENTIFICATION.R
Notes: This R script (Script 12) performs SNP-based species identification by comparing the species-specific SNPs ascertained in Script 9 with the SNPs called from sediment samples in Script 11.
Variant positions, reference and alternative alleles, and read depth counts per sample were extracted from the each genus-specific SNP dataset using R (R Core Team 2022). 
To compare the genus-specific SNPs from sediment samples with ascertained species-specific SNPs from the reference alignment, the two datasets were merged in R (R Core Team 2022) based on genomic position and allele identity. Genus-specific SNPs that matched ascertained species-specific SNPs were assigned to the corresponding species (Appendix S2: Table S7), whereas SNPs absent from all reference genomes were classified as "unknown." Variants occurring at positions where all reference genomes were identical were excluded from further analysis, as they were not phylogenetically informative. The absolute read counts and the relative proportion of reads assigned to each reference species were calculated and visualized for each site (Appendix S1: Figure S20-S29). Data were analyzed and visualized in R (R Core Team 2022) using the package tidyverse (Wickham et al. 2019).

Script: 13.1_MAPDAMAGE_PINACEAE_READS.sl
Note: Script to process all sediment cores simultaneously. It assesses ancient DNA damage patterns from Pinaceae classified reads mapped to the concatenated reference genome, by merging all samples (sorted BAM files) per-site.
Script: 13.2_MAPDAMAGE_GENUS_READS.sl
Note: Script for processing a single sediment core (e.g., Khamra), to be adapted for each individual core. It assesses ancient DNA damage patterns from genus-specific reads mapped to distant species reference genomes, by merging all samples (sorted BAM files) per site, separately for each genus.
Script: 13.3_MAPDAMAGE_PLOT_RESULTS.R
Note: Script to plot the results from all sediment cores simultaneously.
Ancient damage patterns were assessed from the Pinaceae classified reads mapped to the concatenated genome (Appendix S1: Figure S30, upper panels) and from the genus-specific reads mapped to the distant species reference genome (Appendix S1: Figure S30, lower panels). All sample alignment files (BAM) of each site were merged together using SAMTOOLS MERGE (H. Li and Durbin 2009). Ancient damage patterns were assessed using MAPDAMAGE2 (v. 2.2.1) (Jónsson et al. 2013) providing the concatenated genome and the distant species reference genomes used for sample alignments (Appendix S2: Table S3A and 3C). The frequency of C>T substitutions across the first 25 positions at the 5' end of both Pinaceae classified reads, and genus-specific reads were extracted from the 5pCtoT_freq output file and plotted for each site.

Script: 14.1_REMOVE_CT_FROM_VCF.sl
Note: Script to remove C-T transitions from genus-specific VCF files prior to admixture and TCS network analyses.
Script: 14.2_ADMIXTURE_ANALYSIS.sl and 14.2_PLINK_PARAMETERS.txt
Note: Scripts to perform admixture analysis using genus-specific VCF files, to be adapted for each specific genus.
Script:14.3_PLOT_CROSS_VALIDATION_RESULTS.R and 14.3_PLOT_ADMIXTURE_RESULTS.R
Note: Scripts for plotting genus-specific admixture analysis results, to be adapted for each specific genus.
Admixture analysis was conducted using genus-specific SNPs after excluding post-mortem derived C>T substitutions, utilizing ADMIXTURE (v. 1.3.0) (Alexander, Novembre, and Lange 2009). To determine the optimal number of clusters (K), cross-validation analysis was performed (Appendix S1: Figure S31), evaluating up to 10 clusters. The results for the two values of (K) with the lowest cross-validation errors (the optimal (K) and the next best (K)) were plotted (Appendix S1: Figure S32-S35). 

Note: TCS networks are reproducible using the NEXUS files for each genus listed in Appendix S2: Table S8, by selecting the samples corresponding to each TCS network figure. 
Population TCS haplotype networks were assessed between sediment samples using genus-specific SNPs. Post-mortem derived C>T substitutions were removed from the genus-specific variant call files (VCF) using BCFTOOLS (v. 1.9) (H. Li 2011). For each Pinaceae genus, the VCF files were converted into NEXUS format using PGDspider with haploid settings (Lischer and Excoffier 2012). The NEXUS files were then used to construct TCS haplotype networks (Templeton, Crandall, and Sing 1992) via POPART v. 1.7 (Clement et al. 2002) with default settings (Appendix S2: Table S8). Only the SNPs with available data in all selected samples were included. Missing data rates for each sample were considered neutral, neither supporting nor contradicting the connections between nodes. Sample selection for each TCS network was guided by the research focus on specific species and their potential links across sites, while also accounting for the fact that including low-coverage samples would significantly reduce the number of informative segregating sites.

References

Alexander, D. H., J. Novembre, and K. Lange. 2009. Fast model-based estimation of ancestry in unrelated individuals. Genome Research 19:1655–1664.

Andrews, S. 2015. FastQC: A quality control tool for high throughput sequence data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

Barnett, D. W., E. K. Garrison, A. R. Quinlan, M. P. Str̈mberg, and G. T. Marth. 2011. Bamtools: A C++ API and toolkit for analyzing and managing BAM files. Bioinformatics 27:1691–1692.

Bushnell, B. 2014. BBMap: A Fast, Accurate, Splice-Aware Aligner.

Chang, C. C., C. C. Chow, L. C. A. M. Tellier, S. Vattikuti, S. M. Purcell, and J. J. Lee. 2015. Second-generation PLINK : rising to the challenge of larger and richer datasets. GigaScience 25:4–7.

Chen, S., Y. Zhou, Y. Chen, and J. Gu. 2018. Fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34:i884–i890.

Clement, M., Q. Snell, P. Walke, D. Posada, and K. Crandall. 2002. TCS: Estimating gene genealogies. Proceedings 16th International Parallel and Distributed Processing Symposium.

Danecek, P., A. Auton, G. Abecasis, C. A. Albers, E. Banks, M. A. Depristo, R. E. Handsaker, G. Lunter, G. T. Marth, S. T. Sherry, G. Mcvean, R. Durbin, G. Project, and T. Vcf. 2011. The variant call format and VCFtools 27:2156–2158.

Garrison, E., and G. Marth. 2012. Haplotype-based variant detection from short-read sequencing.

Li, H., B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, and R. Durbin. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25:2078–2079.

Li, H., and R. Durbin. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754–1760.

Lischer, H. E. L., and L. Excoffier. 2012. PGDSpider: An automated data conversion tool for connecting population genetics and genomics programs. Bioinformatics 28:298–299.

Lockwood, J. D., J. M. Aleksić, J. Zou, J. Wang, J. Liu, and S. S. Renner. 2013. A new phylogeny for the genus Picea from plastid, mitochondrial, and nuclear sequences. Molecular Phylogenetics and Evolution 69:717–727.

Lu, J., N. Rincon, D. E. Wood, F. P. Breitwieser, C. Pockrandt, B. Langmead, S. L. Salzberg, and M. Steinegger. 2022. Metagenome analysis using the Kraken software suite. Nature Protocols 17:2815–2839.

O’Leary, N. A., M. W. Wright, J. R. Brister, S. Ciufo, D. Haddad, R. McVeigh, B. Rajput, B. Robbertse, B. Smith-White, D. Ako-Adjei, A. Astashyn, A. Badretdin, Y. Bao, O. Blinkova, V. Brover, V. Chetvernin, J. Choi, E. Cox, O. Ermolaeva, C. M. Farrell, T. Goldfarb, T. Gupta, D. Haft, E. Hatcher, W. Hlavina, V. S. Joardar, V. K. Kodali, W. Li, D. Maglott, P. Masterson, K. M. McGarvey, M. R. Murphy, K. O’Neill, S. Pujar, S. H. Rangwala, D. Rausch, L. D. Riddick, C. Schoch, A. Shkeda, S. S. Storz, H. Sun, F. Thibaud-Nissen, I. Tolstoy, R. E. Tully, A. R. Vatsan, C. Wallin, D. Webb, W. Wu, M. J. Landrum, A. Kimchi, T. Tatusova, M. DiCuccio, P. Kitts, T. D. Murphy, and K. D. Pruitt. 2016. Reference sequence (RefSeq) database at NCBI: Current status, taxonomic expansion, and functional annotation. Nucleic Acids Research 44:D733–D745.

Oliva, A., R. Tobler, A. Cooper, B. Llamas, and Y. Souilmi. 2021. Systematic benchmark of ancient DNA read mapping. Briefings in Bioinformatics 22:1–12.

Templeton, A. R., K. A. Crandall, and C. F. Sing. 1992. A cladistic analysis of phenotypic associations with haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation. Genetics 132:619–633.

Wickham, H., M. Averick, J. Bryan, W. Chang, L. D. A. Mcgowan, R. François, G. Grolemund, A. Hayes, L. Henry, J. Hester, M. Kuhn, T. Lin, E. Miller, S. M. Bache, K. Müller, J. Ooms, D. Robinson, D. P. Seidel, V. Spinu, K. Takahashi, D. Vaughan, C. Wilke, and K. Woo. 2019. Welcome to the Tidyverse Tidyverse. The Journal of Open Source Software 4:1686.

Wood, D. E., J. Lu, and B. Langmead. 2019. Improved metagenomic analysis with Kraken 2. Genome Biology 20:257.
