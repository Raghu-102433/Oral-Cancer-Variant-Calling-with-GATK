🧬 Oral Cancer Variant Calling Using WES

This project identifies and interprets germline and somatic variants from whole exome sequencing (WES) data of oral cancer samples. The workflow integrates data preprocessing, variant calling with GATK, functional annotation using Ensembl VEP, and comparative analysis.

Two matched Illumina samples were used:

SRR32583844 — Healthy (normal/control)

SRR32583845 — Tumour (oral squamous cell carcinoma)

⚙️ Getting Started
```bash
Step 1: Create Conda Environment

conda env create -f environment.yml

conda activate ngs_env

Step 2: Run Pipelines

bash germline_variantcalling_pipeline.sh

bash somatic_variantcalling_pipeline.sh
```

📊 Summary

The processed data was analyzed to identify variants with high functional impact:

Germline analysis revealed distinct sets of high-impact variants between healthy (SRR32583844) and tumor (SRR32583845) samples. Some variants were shared but showed increased frequency in tumor tissue, indicating possible clonal expansion.

Somatic analysis uncovered mutations enriched in immune checkpoint genes (e.g., LILRB1/2, HLA loci), mucin genes (MUC3A, MUC12), and chromatin remodelers (KDM6A, KMT2C). These support known hallmarks of oral cancer including immune evasion and epigenetic dysregulation.

These findings contribute to understanding variant patterns in oral cancer genomes and highlight candidate genes for future investigation or screening.

📬 Contact
For questions, suggestions, or collaborations, please contact:
📧 raghava.332410@gmail.com

© 2025 | WES-Based Oral Cancer Genomics Project
