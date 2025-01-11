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
Colocalization analysis was conducted using the Coloc R package (version 5.2.3) in R (version 4.4.2). The analysis was executed with the following parameters for the SNP priors:

p1 = 1e-04
p2 = 1e-04
p12 = 1e-05

```
gwas <- read.table("rs1385742_gwas.txt", header=TRUE)
eqtl <- read.table("CD2AP_eqtl.txt", header=TRUE)
rownames(gwas) <- gwas$SNP; rownames(eqtl) <- eqtl$SNP

my.res <- coloc.abf(
  dataset1 = list(snp = gwas$SNP, pvalues = gwas$Pvalue, beta = gwas$Beta, varbeta = gwas$Varbeta, MFA = gwas$MFA, position = gwas$Location, N = gwas$Size, type = "cc"),
  dataset2 = list(snp = eqtl$SNP, pvalues = eqtl$Pvalue, beta = eqtl$Beta, varbeta = eqtl$Varbeta, MFA = eqtl$MFA, position = eqtl$Location, N = eqtl$Size, type = "cc")
)

# Sorting the results by posterior probability of H4 (colocalization)
o <- order(my.res$results$SNP.PP.H4, decreasing=TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp

results <- my.res$results
write.table(results, file = "rs1385742_CD2AP_SCVA_01062025.txt", sep = "\t", quote = FALSE)
```
# Results
The Approximate Bayes Factor (ABF) colocalization analysis revealed strong evidence for colocalization in the regions of rs4311 and rs1385742. The posterior probability for H4 (indicating a shared causal variant) was high for both of these regions.
| **SNP**    | **H0**    | **H1**    | **H2**   | **H3**   | **H4**   |
| ---------- | --------- | --------- | -------- | -------- | -------- |
| rs4311     | 9.83E-68  | 1.35E-65  | 2.03E-04 | 2.69E-02 | 9.73E-01 |
| rs1385742  | 1.60E-18  | 2.32E-13  | 6.80E-08 | 8.85E-03 | 9.91E-01 |
| rs1859788  | 2.37E-295 | 1.11E-283 | 2.13E-12 | 1.00E+00 | 1.33E-14 |
| rs12151021 | 4.74E-31  | 2.31E-24  | 2.06E-07 | 1.00E+00 | 9.58E-06 |

The analysis also identified the specific SNPs that are likely associated with both traits in the rs4311 and rs1385742 regions. The highest posterior probability was found for rs1385742 as the shared causal variant.
| **snp**   | **position** | **SNP.PP.H4** |
| --------- | ------------ | ------------- |
| rs4291    | 63476833     | 0.9418704     |
| rs4292    | 63476980     | 0.049837      |
| rs1385742 | 47627419     | 0.998549155   |

Table: Coloc based colocalization analysis results of lead variants.
| **snp**    | **position** | **SNP.PP.H4** |
| ---------- | ------------ | ------------- |
| rs4311     | 63483402     | 2.53E-33      |
| rs1385742  | 47627419     | 0.998549155   |
| rs1859788  | 100374211    | 4.27E-257     |
| rs12151021 | 1050875      | 2.92E-21      |

### END ###






