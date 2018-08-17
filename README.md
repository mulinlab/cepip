# cepip
Context-dependent epigenomic weighting for regulatory variant prioritization


# Introduction

Majority of trait/disease associated variants identified by genome wide association studies (GWASs) locate in the regulatory regions. Since gene regulation is highly context-specific, it remains challenging to fine-map and prioritize functional regulatory variants in a particular cell/tissue type and apply them to disease-associated genes detection. By connecting large-scale epigenome profiles to expression quantitative trait loci (eQTLs) in a wide range of human tissues/cell types, we identify combination of several critical chromatin features that predict variant regulatory potential. We develop a joint likelihood framework to measure the regulatory probability of genetic variants in a context-dependent manner. We show our model is superior to existing cell type-specific methods and exhibit significant GWAS signal enrichment. Using phenotypically relevant epigenomes to weight GWAS SNPs, we discover more disease-associated genes owing to regulatory changes and improve the statistical power in gene-based association test.

# Installing cepip

Simply decompress the archive and run the following command.

   java -Xms256m  -Xmx1300m -jar ./cepip.jar  [arguments] 

The arguments -Xms256m and -Xmx1300m set the initial and maximum Java heap sizes for cepip as 256 megabytes and 1.3 gigabytes respectively. Specifying a larger maximum heap size can speed up the analysis. A higher setting like -Xmx2g or even-Xmx5g is required when there is a large number of variants, say 5 million. The number, however, should be less than the size of physical memory of a machine.

Note: [arguments] can be saved in a flat text file.

More infomation about installation, usage and updates, please visit: http://147.8.193.36/cepip
