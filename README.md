# Evolution and plasticity of gene expression under progressive warming in Drosophila subobscura
Authors: Marta A. Antunes, Marta A. Santos, Ana S. Quina, Mauro Santos, Margarida Matos & Pedro Simões
# - Workflow -

Below we provide a step-by-step workflow of the analyses in the paper. Shaded background represents the code for the analysis or filenames for specific files containing code (provided).  

DESeq2, a widely used tool for differential gene expression analysis, estimates variance-mean dependence in count data from high-throughput sequencing assays and tests for differential expression using the negative binomial distribution. However, DESeq2 does not provide p-values for interaction terms in its analyses.  

To address this limitation, some of the code in this repository (e.g. run_Overall_gene_expression_analysis.R) extends DESeq2's functionality. It enables differential gene expression analysis based on the negative binomial distribution while incorporating interaction terms, allowing for more nuanced statistical insights.

# 1. Pre-processing
1.1. Quality control (fastqc is a quality control tool for high throughput sequence data)  
input files: the “.fq.gz” files stored at Sequence Read Archive (SRA) with accession number: XXX 
```
for file in  /*.fq.gz; do
fastqc $file -t 11 -q -o /1_preprocessing/; done
```
output files: “.fq_fastqc” file for each “.fq” file 

1.2. Trimming  
input files: the “.fq.gz” files stored at Sequence Read Archive (SRA) with accession number: XXX 
```
fastp -i ${file}_1.fq.gz -I ${file}_2.fq.gz -o /2_preprocessing/${file}_1.Q20L120.fq.gz -O 2_preprocessing/${file}_2.Q20L120.fq.gz --unpaired1 /2_preprocessing/${file}_unpaired1 --unpaired2 /2_preprocessing/${file}_unpaired2 --thread 11 --adapter_sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --average_qual 20 --length_required 120 --overrepresentation_analysis  --json /1_preprocessing/${file}.json --html /1_preprocessing/${file}.html
```
outputs: “.Q20L120.fq.gz” file for each input file

# 2. Mapping
2.1 Index genome and annotation file
```
/opt/tools/STAR-2.7.9a/source/STAR --runMode genomeGenerate --runThreadN 11 --genomeDir /3_mapping/StarIndexFiles/ --genomeFastaFiles /ReferenceGenome/drosophilaSubobscura_1.0/ncbi/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna --sjdbGTFfile /ReferenceGenome/drosophilaSubobscura_1.0/ncbi/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.gff --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene Parent --genomeSAindexNbases 12
```

2.2 Pass 1 for one sample (as example), will generate a Sam file and one SJ.tab file for each pair
```
/opt/tools/STAR-2.7.9a/source/STAR --runMode alignReads --runThreadN 11 --genomeDir /3_mapping/StarIndexFiles/ --readFilesCommand zcat --readFilesIn PT1_G23_C_1.Q20L120.fq.gz PT1_G23_C_2.Q20L120.fq.gz --outFileNamePrefix RNASEQ_PT1G23C.Q20L120_Alingment1
```

2.3 Pass 2 for one sample (as example), give as argument the SJ.tab file produced for all the pair. Generates a sorted by coordinate BAM file
```
/opt/tools/STAR-2.7.9a/source/STAR --runMode alignReads --runThreadN 11 --genomeDir /3_mapping/StarIndexFiles/ --readFilesCommand zcat --readFilesIn PT1_G23_C_1.Q20L120.fq.gz PT1_G23_C_2.Q20L120.fq.gz --sjdbFileChrStartEnd *SJ.out.tab --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNASEQ_PT1G23C.Q20L120_Alingment2 --limitBAMsortRAM 3059132604
```

# 3. Feature counts
3.1 run featureCounts  
input files: “.bam” files from the Pass 2 alignment and gff from the reference genome
```
for file in folder
do;/opt/tools/subread-2.0.0-Linux-x86_64/bin/featureCounts -p -T 9 -F GFF -t gene -g ID -s 2 -C -a  /ReferenceGenome/drosophilaSubobscura_1.0/ncbi/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.gff -o /4_featureCounts/${file}_counts  /3_mapping/STAR/2ndRun/RNASEQ_${file}.Q20L120_Alingment2Aligned.sortedByCoord.out.bam > ${file}.log 2>&1;done
```
output files: count files for each sample

3.2 generate supplementary figure S2  
MultiQC software receives count files for all samples and generates the figure. The dot below (.) indicates that the software was run in the current working directory.
```
Multiqc .

```

# 4. Overall gene expression analysis
4.1 normalize counts with DESeq2 within galaxy  
input files: count files for all samples  
output file: Galaxy238.tabular (file with normalized counts)  
all normalizations were done selection the option “output options”>”Output normalised counts” on DESeq2 within Galaxy

4.2 make PCA and PVCA for all samples using our full gene dataset  
input file: Galaxy238.tabular
```
normalizedCounts2PCA.R
```
output: PCA and PVCA plots (Fig. 2A)

4.3 differential expression analysis  
input file: Galaxy238.tabular 
```
run_Overall_gene_expression_analysis.R
```
output file: results_OverallGeneExpressionAnalysis.csv (supplementary tables S9 and S10)

# 5. Evolutionary changes in the warming environment 
5.1 normalize counts with DESeq2 within galaxy for each latitudinal population
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
| Input file                     | count files for PT_W and WPT_W samples | count files for NL_W and WNL_W samples |
| Output file                    | /Analysis_of_Selection/Galaxy218-Normalized_counts.tabular | /Analysis_of_Selection/Galaxy225-Normalized_counts.tabular |

5.2 differential expression analysis in each latitudinal population

`/Analysis_of_Selection/run_Analysis_of_Selection.R` (select appropriate file names within the code)

|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
| Input file	                 |/Analysis_of_Selection/Galaxy218-Normalized_counts.tabular|	/Analysis_of_Selection/Galaxy225-Normalized_counts.tabular|
|Output file	                 |supplementary table S12A	   |supplementary table S12B         |

5.3 make Manhattan plot for each latitudinal population

`/Manhattan plots/make_plot_with_pvalues_negativelog2FC_below.R` (select appropriate file names within the code)

|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input file	                 |/Manhattan plots/inputFilePT.csv |	/Manhattan plots/inputFileNL.csv|
|Output file	                 |Fig. 3A	                   |Fig. 3B

5.4 make PCA and PCVA for both latitudinal populations when tested in the warming environment using our full gene dataset  
input file: Galaxy238.tabular
```
normalizedCounts2PCA_onlyWenv.R
```
output: PCA and PVCA plots in Fig. 2B

5.5 make PCA and PCVA for the set of candidate genes of each latitudinal population  
`/PCAs/normalizedCounts2PCA_onlyWenv_candidateGenes.R` (select appropriate file names within the code)

|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files                     |	Galaxy238.tabular+/PCAs/CandidatesLowLat.txt|	Galaxy238.tabular+/PCAs/CandidatesHighLat.txt|
|Output	                         |PCA and PVCA Fig. S6A            |	PCA and PVCA Fig. S6B        |


5.6 make Gene Ontology analysis

`/Analysis_of_Selection /prepare_GO_analysis.py` (select appropriate file names within the code)

|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |GCF_008121235.1_UCBerk_Dsub_1.0_genomic.gff+ /Analysis_of_Selection/CandidatesLowLat_noprefix.txt |	GCF_008121235.1_UCBerk_Dsub_1.0_genomic.gff+ /Analysis_of_Selection/CandidatesHighLat_ noprefix.txt|



this python script outputs D. subobscura protein ids for candidate genes in each latitudinal population
detect D. melanogaster orthologs for all D. subobscura proteins (with Proteinortho within galaxy)
input for Proteinortho: GCF_008121235.1_UCBerk_Dsub_1.0_protein.faa and GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa
output from Proteinortho: Galaxy272.tabular D. melanogaster proteins
after excluding non-candidate D. melanogaster protein IDs, the remaining relevant IDs were used as input for batch Entrez.
the gene ids from batch entrez were used in BiNGO app from Cytoscape to make the enrichment analysis. The parameters used were the default ones with two exceptions: Ontology file that was set to GO_Biological_Process, GO_Molecular_Function or GO_Cellular_Component and organism/annotation was set to D. melanogaster.

5.7 Analysis of the effects of selection in the magnitude of differences between latitudinal populations  
5.7.1 normalize counts with DESeq2 within galaxy  
input files: count files for PT_W, NL_W, WPT_W and WNL_W samples  
output file: /Analysis_of_Convergence/Galaxy232-Normalized_counts.tabular

5.7.2 differential expression analysis  
“input_analysis4.csv” file was created by editing the “Galaxy232-Normalized_counts.tabular” file

`/Analysis_of_Convergence/run_Analysis_convergence.R` (select appropriate file names within the code)


|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |/Analysis_of_Convergence/input_analysis4.csv+ /Analysis_of_Convergence/CandidatesLowLat.txt|	/Analysis_of_Convergence/input_analysis4.csv+ /Analysis_of_Convergence/CandidatesHighLat.txt|
|Output                          |	Table S14A                 |	Table S14B                   |



# 6. Evolution of plasticity
6.1 Plasticity of the control populations  
6.1.1 normalize counts with DESeq2 within galaxy for each latitudinal population
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |count files for PT_C and PT_W samples	|count files for NL_C and NL_W samples|
|Output                          |	/Analysis_of_Plasticity/Galaxy108-Normalized_counts.tabular	|/Analysis_of_Plasticity/Galaxy123-Normalized_counts.tabular|

6.1.2 differential expression analysis on each latitudinal population  
`/Analysis_of_Plasticity/run_Analysis_of_Plasticity.R` (select appropriate file names within the code)
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |/Analysis_of_Plasticity/Galaxy108-Normalized_counts.tabular	|/Analysis_of_Plasticity/Galaxy123-Normalized_counts.tabular|
|Output                          | Supplementary Table S15A	   |Supplementary Table S15B         |



6.2 Analysis of the effects of history on the plasticity of the control populations  
6.2.1 normalize counts with DESeq2 within galaxy  
input file: count files for PT_C, NL_C, PT_W and NL_W populations  
output file: /EffectOfHistoy_inPlasticityOfControlPops/Galaxy179-Normalized_counts.tabular

6.2.2 differential expression analysis  
input file: /EffectOfHistoy_inPlasticityOfControlPops/Galaxy179-Normalized_counts.tabular
```
/EffectOfHistoy_inPlasticityOfControlPops/run_Analysis_EffectOfHistoryInPlastOFControl.R
```
output file: Supplementary Tables S16A and S16B

6.3 Analysis of differences between control populations of distinct history  
6.3.1 normalize counts with DESeq2 within galaxy for each environment
|                                | Control environment        | Warming environment       |
|--------------------------------|----------------------------|---------------------------|
|Input files	                 |count files for PT_C and NL_C samples|	count files for PT_W and NL_W samples|
|Output	                         |/Analysis_of_History/Galaxy377-Normalized_counts.tabular	|/Analysis_of_History/Galaxy380-Normalized_counts.tabular|

6.3.2 Differential expression analysis for each environment  
`/Analysis_of_History/run_Analysis_of_History.R` (select appropriate file names within the code)
|                                | Control environment        | Warming environment       |
|--------------------------------|----------------------------|---------------------------|
|Input files	                 |/Analysis_of_History/Galaxy377-Normalized_counts.tabular|/Analysis_of_History/Galaxy380-Normalized_counts.tabular|
|Output                          |Supplementary Table S17A    |	Supplementary Table S17B  |



6.4 Plasticity of the Warming populations
6.4.1 normalize counts with DESeq2 within galaxy for each latitudinal population
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |count files for WPT_C and WPT_W samples|	count files for WNL_C and WNL_W samples|
|Output                  	 |/Analysis_of_Plasticity/Galaxy359-Normalized_counts.tabular	|/Analysis_of_Plasticity/Galaxy362-Normalized_counts.tabular|

6.4.2 Differential expression analysis
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |/Analysis_of_Plasticity/Galaxy359-Normalized_counts.tabular|	/Analysis_of_Plasticity/Galaxy362-Normalized_counts.tabular|
|Output	                         |Supplementary Table S18A	   |Supplementary Table S18B         |

```
/Analysis_of_Plasticity/run_Analysis_of_Plasticity.R
```
 (select appropriate file names within the code)

6.5 Evolution of Plasticity
6.5.1 normalize counts with DESeq2 within galaxy for each latitudinal population
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |count files for PT_C, WPT_C, PT_W and WPT_W samples|	count files for NL_C, WNL_C, NL_W and WNL_W samples|
|Output                    	 |/Analysis_of_the_Evolution_of_plasticity/Galaxy255-Normalized_counts.tabular|	/Analysis_of_the_Evolution_of_plasticity/Galaxy252-Normalized_counts.tabular|

6.5.2 Differential expression analysis
	Low latitude populations	High latitude populations
Input files	/Analysis_of_the_Evolution_of_plasticity/Galaxy255-Normalized_counts.tabular	/Analysis_of_the_Evolution_of_plasticity/Galaxy252-Normalized_counts.tabular
Output	Supplementary Table S19A	Supplementary Table S19B

```
/Analysis_of_the_Evolution_of_plasticity/run_Analysis_of_evolution_of_plasticity.R
```
 (select appropriate file names within the code)

# 7. Plasticity vs Evolution (files provided)
7.1 Fisher exact test for each latitudinal population
input file: not necessary
```
Fisher_test_for_plasticity_and_selection.R
```
output: p-value, alternative hypothesis, 95 percent confidence interval and odds ratio for each latitudinal population

# 8. Jaccard index (files provided)
|                                | Low latitude populations        | High latitude populations       |
|--------------------------------|---------------------------------|---------------------------------|
|Input files	                 |/Analysis_of_Selection/Galaxy218-Normalized_counts.tabular+interest.txt|	/Analysis_of_Selection/Galaxy225-Normalized_counts.tabular+interest.txt|
|Output	                         |Mean jaccard Index	           |Mean jaccard Index               |

```
calculate_jaccardIndex.R
```

