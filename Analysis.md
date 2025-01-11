# Plan of Action

1. **Download the GWAS and eQTL datasets** from the respective websites.
2. **LiftOver the GWAS dataset** from hg19 to hg38 genomic coordinates.
3. **Subset the GWAS data** to extract a region centered around the lead SNPs.
4. **Subset the eQTL data** for genes associated with the lead variants.
5. **Run the Coloc analysis** and interpret the results.

# Downloading the GWAS Data

The GWAS summary statistics were downloaded from the UK EBI FTP server using the `wget` command:

```
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012877/GCST90012877_buildGRCh37.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012877/GCST90012877_buildGRCh37.tsv.gz-meta.yaml
```

# Liftover of GWAS Data (hg19 to hg38)
To convert the genomic coordinates of the GWAS data from hg19 to hg38, the UCSC LiftOver tool was used. The following steps were followed:

## Download the hg19ToHg38.over.chain file:
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```
## Convert the GWAS data to BED format using awk:
```
awk 'BEGIN {OFS="\t"} NR > 1 {print "chr"$3, $4-1, $4, $0}' GCST90012877_buildGRCh37.tsv > lift_input.bed
```
## Perform the Liftover conversion:
```
liftOver -bedPlus=3 -tab lift_input.bed hg19ToHg38.over.chain lift_output.bed lift_unmapped.bed
```
## Reformat the output and merge the converted coordinates into the GWAS dataset:
```
awk 'BEGIN {OFS="\t"} {$7 = $3; $6 = substr($1, 4); print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23}' lift_output.bed > GCST90012877_buildGRCh38.tsv
```
Note: A total of 4,972 SNPs were lost during the Liftover process, likely due to sequence differences between assemblies. However, the presence of lead variants was manually verified.

# Downloading eQTL Data
Since the lead variants were located on chromosomes 6, 7, 17, and 19, eQTL data from these chromosomes were sufficient for the analysis. The relevant eQTL datasets were downloaded using the following wget commands:
```
wget https://download.metabrain.nl/2021-07-23-release/2021-07-23-cortex-EUR-80PCs-chr6.txt.gz
wget https://download.metabrain.nl/2021-07-23-release/2021-07-23-cortex-EUR-80PCs-chr7.txt.gz
wget https://download.metabrain.nl/2021-07-23-release/2021-07-23-cortex-EUR-80PCs-chr17.txt.gz
wget https://download.metabrain.nl/2021-07-23-release/2021-07-23-cortex-EUR-80PCs-chr19.txt.gz
```

# Extracting GWAS and eQTL Data for Analysis
The coordinates and information for each lead SNP were extracted from the Liftover GWAS dataset using a custom Perl script. SNPs within a 500KB window upstream and downstream of each lead SNP were then selected.

The following columns were extracted from the GWAS summary statistics:

SNP: Variant ID
CHR: Chromosome
Location: Base pair location
Pvalue: P-value
MFA: Effect allele frequency
Beta: Beta coefficient
SE: Standard error
Sample Size: 472,868
Similarly, for the eQTL dataset, the following columns were extracted:

SNP: SNP identifier
MetaP: Meta P-value
MFA: SNP effect allele frequency
MetaBeta: Meta beta coefficient
MetaSE: Meta standard error
MetaPN: Meta sample size
The data were stored in tab-delimited files for each gene and lead variant. An example of the Perl command used to extract data:

```
perl ProcessGWAS_EQTL_V3.pl rs1859788 GCST90012877_buildGRCh38.tsv PILRA 2021-07-23-cortex-EUR-80PCs-chr7.txt
```


# Colocalization Analysis
