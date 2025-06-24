# 🧬 Somatic Variant Analysis Summary (Healty vs Tumour)

This summary describes somatic mutations identified from the tumor sample SRR32583845 using GATK Mutect2, contrasted against the matched normal sample SRR32583844.

## Summary Statistics

- Total raw somatic variants: 78,214
- Filtered high-confidence somatic variants: 12,603
- High-impact variants: 1,684 (stop-gain, frameshift, etc.)

## Top Somatically Mutated Genes (High Impact)

From the filtered high-impact somatic variants, the following genes were most frequently affected:

- FIP1L1 (n = 36): oncogenic fusion gene known in leukemias
- LILRB1 & LILRB2 (n = 29, 15): immune checkpoint regulators
- HLA-A/B/C/D loci: immunogenic variation; reflects immune escape
- MUC3A, MUC6, MUC12: mucin genes altered in epithelial cancers
- KDM6A, KMT2C: chromatin remodeling genes frequently mutated in tumors

Additional frequent variants were found in pseudogenes and immune-related regulatory regions.

##  Interpretation

The mutation landscape reveals:

- Extensive changes in immune-related loci, particularly in HLA and LILRB families, suggesting immune evasion mechanisms.
- Dysregulation of mucins and epithelial defense genes (e.g., MUC3A, MUC12) supports a mucosa-origin tumor microenvironment.
- Mutations in chromatin modifiers (KDM6A, KMT2C) indicate possible epigenetic dysregulation.

These patterns reflect a tumor profile characteristic of epithelial malignancies with a strong immunomodulatory component.

## Conclusion

Somatic variant analysis of SRR32583845 (tumor) against SRR32583844 (normal) revealed a diverse panel of high-impact mutations affecting immune checkpoints, chromatin remodelers, and epithelial defense genes. These alterations may underpin immune escape, genomic instability, and selective clonal expansion—hallmarks of oral cancer pathogenesis.
