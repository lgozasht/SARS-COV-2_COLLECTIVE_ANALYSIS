# SARS-COV-2_COLLETIVE_ANALYSIS

We previously developed a program in python3 that flags possible lab-associated variants in SARS-CoV-2 data (see https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter). This repository hosts our improved pipeline which employs similar methods as *COVID-19-Lab-Specific-Bias-Filter* but includes additional modifications. 

We describe these modifications in detail in https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/12.

Briefly, our improved method requires an "*unresolved*" VCF as input, in which ambiguities remain unresolved (ambiguities can be a valuable signature of potentially errant alleles). It then performs two independent associations: one after resolving all ambiguous calls as the alternative allele (“alternative resolved”), and one where all ambiguous calls are resolved by parsimony (selecting the reference allele where there are equally parsimonious configurations: “parsimony resolved”). In addition, we added a “country-association” module, using the same methods we previously employed for detecting “lab-associated” variants. Our pipeline also associates **specific** alleles (including ambiguities) with specific labs and countries (i.e. if 99% of "W" allele calls at a given site are attributed to a particular lab, regardless if that lab contributes significantly to the total alternate allele count).

##Usage
