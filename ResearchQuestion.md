# Research Question

The goal of this analysis is to perform a **Coloc**-based colocalization analysis using the **Coloc** R package to identify relationships between the genetic regulation of expression levels (eQTLs) and local Alzheimer’s disease (AD) GWAS signals.

## Specific Instructions

1. Perform colocalization analysis using the **Coloc** R package.
   - Use a 1-Megabase window range centered on the GWAS lead variants.
   - Assume a single causal variant for each region.
   - For the eQTL data, use only the "Cortex-EUR" data from the listed resource.

2. Evaluate the following lead variants and their nearby genes:
   - **rs4311** (ACE)
   - **rs1385742** (CD2AP)
   - **rs1859788** (PILRA)
   - **rs12151021** (ABCA7)

## Resources

The necessary datasets and tools are publicly available:

- **Coloc R package**: [https://cran.r-project.org/web/packages/coloc/index.html](https://cran.r-project.org/web/packages/coloc/index.html)
- **GWAS Summary Statistics (hg19)**: [https://www.nature.com/articles/s41588-020-00776-w](https://www.nature.com/articles/s41588-020-00776-w)
  - Data Identifier: GCST90012877
- **Brain eQTL Summary Statistics (hg38)**: [https://www.nature.com/articles/s41588-023-01300-6](https://www.nature.com/articles/s41588-023-01300-6)
- **Metabrain eQTL Data**: [https://download.metabrain.nl/files.html](https://download.metabrain.nl/files.html)

