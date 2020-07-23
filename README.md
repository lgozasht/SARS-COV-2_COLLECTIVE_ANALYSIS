# SARS-COV-2_COLLECTIVE_ANALYSIS

We previously developed a program in python3 that flags possible lab-associated variants in SARS-CoV-2 data (see https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter). This repository hosts our improved pipeline which employs similar methods as *COVID-19-Lab-Specific-Bias-Filter* but includes additional modifications. 

We describe these modifications in detail in https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12.

Briefly, our improved method requires an "*unresolved*" VCF as input, in which ambiguities remain unresolved (ambiguities can be a valuable signature of potentially errant alleles). It then performs two independent associations: one after resolving all ambiguous calls as the alternative allele (“alternative resolved”), and one where all ambiguous calls are resolved by parsimony (selecting the reference allele where there are equally parsimonious configurations: “parsimony resolved”). In addition, we added a “country-association” module, using the same methods we previously employed for detecting “lab-associated” variants. Our pipeline also associates **specific** alleles (including ambiguities) with specific labs and countries (i.e. if 99% of "W" allele calls at a given site are attributed to a particular lab, regardless if that lab contributes significantly to the total alternate allele count). 

We also provide a separate program that identifies local positive linkage disequilibrium (LD) between variants. Local positive LD can be a signature of recurrent error.

## Dependencies

#### strain_phylogenetics 

Install from https://github.com/yatisht/strain_phylogenetics

#### PLINK 2

Install from https://www.cog-genomics.org/plink/2.0/ and add to PATH

## Association Usage

SARS-COV-2_COLLETIVE_ANALYSIS.py only requires a GISAID metadata file, *unresolved* VCF file, and corresponding newick tree as input. The pipeline will automatically "filter" the input VCF file for samples that also exist in the corresponding tree and will produce "filtered_unresolved.vcf" as a result. It will then produce "filtered_resolved_alt.vcf" after resolviong ambiguities as alternate alleles and "filtered_resolved_ref.vcf" after resolving ambiguities using parsimony. Note that the pipeline will not overwrite these files, so to rerun on a new dataset just delete them or remove them from your current working directory.

```
python3 SARS-COV-2_COLLETIVE_ANALYSIS.py [options] -m [Path to GISAID metadata file] -v [Path to VCF file] -tree [Path to newick tree] -o [Path to output directory]
```

### Input:

GISAID metadata file, "*unresolved*" VCF and coorresponding newick tree

### Options:

**-min_parsimony *I***: Minimum parsimony (must be an integer) default = 4

**-dependencies**: Check for dependencies

### Output:

**final_table.tsv**:

| Column | Description |
| ------ | ----------- |
| Reference | reference genome  |
| Start | site start |
| Stop | site stop |
| Snp | snp |
| Alt Resolved Parsimony | "alternative resolved" parsimony score |
| MAC | "alternative resolved" minor allele count |
| MAF | "alternative resolved" minor allele frequency |
| Primer Overlap | ARTIC primer overlapping site (NA if none)  |
| Primer Vicinity | ARTIC primer within 10bp of site (NA if none) |
| Country | "alternative resolved" country associated with variant (NA if none) |
| Country Association | "alternative resolved" percent of minor alleles attributed to that country |
| Lab | "alternative resolved" lab associated with variant (NA if none) |
| Lab Association | "alternative resolved" percent of minor alleles attributed to that lab |
| Ref Resolved Parsimony | "parsimony resolved" parsimony score |
| MAC | "parsimony resolved" minor allele count |
| MAF | "parsimony resolved" minor allele frequency |
| Country | "parsimony resolved" country associated with variant (NA if none) |
| Country Association | "parsimony resolved" percent of minor alleles attributed to that country |
| Lab | "parsimony resolved" lab associated with variant (NA if none) |
| Lab Association | "parsimony resolved" percent of minor alleles attributed to that lab |
| Specific Alt Associations | Specific alt associations populate the remaining columns (depending on how many alternate alleles exist at a given site) |


**samples_in_latest_tree.txt**: List of samples that exist both in the provided VCF and the corresponding newick tree.

**filtered_unresolved.vcf**: Unresolved VCF containing only samples that exist both in the provided VCF and the corresponding newick tree.

**filtered_unresolved.txt**: Parsimony file corresponding to the unresolved VCF.

**filtered_resolved_alt.vcf**: "Alternate resolved" VCF

**filtered_resolved_alt.txt**: "Alternate resolved" parsimony file

**filtered_resolved_ref.vcf**: "Parsimony resolved" VCF

**filtered_resolved_ref.txt**: "Parsimony resolved" parsimony file

## Linkage Disequilibrium Usage

Calculates local R^2 to identify linked variants.

```
python3 local_LD.py [options] -v [Path to VCF file] -table [Path to "final_table.tsv" (produced by SARS-COV-2_COLLETIVE_ANALYSIS.py)]  -o [Path to output directory]
```

### Input:

VCF and "final_table.tsv" (produced by SARS-COV-2_COLLETIVE_ANALYSIS.py)

### Options:

**-range *I***: Max range in which LD will be estimated for a given site (Default = 10bp)

**-dependencies**: Check for dependencies

### Output:

**linked_sites.tsv**: 

| Column | Description |
| ------ | ----------- |
| snp | snp  |
| pos | site position |
| linked site\|r^2 | linked site \| corresponding r^2 value |

### References

* Nicola De Maio, Landen Gozashti, Yatish Turakhia, Conor Walker, Robert Lanfear, Russell Corbett-Detig, and Nick Goldman, Issues with SARS-Cov-2 sequencing data: Updated analysis with data from 12th June 2020, Virological [post](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12) 2020.

* Yatish Turakhia, Bryan Thornlow, Landen Gozashti, Angie S. Hinrichs, Jason D. Fernandes, David Haussler, and Russell Corbett-Detig, "Stability of SARS-CoV-2 Phylogenies", bioRxiv [pre-print](https://www.biorxiv.org/content/10.1101/2020.06.08.141127v1) 2020.

