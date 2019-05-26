# cepip
Context-dependent epigenomic weighting for regulatory variant prioritization


# Introduction

Majority of trait/disease associated variants identified by genome wide association studies (GWASs) locate in the regulatory regions. Since gene regulation is highly context-specific, it remains challenging to fine-map and prioritize functional regulatory variants in a particular cell/tissue type and apply them to disease-associated genes detection. By connecting large-scale epigenome profiles to expression quantitative trait loci (eQTLs) in a wide range of human tissues/cell types, we identify combination of several critical chromatin features that predict variant regulatory potential. We develop a joint likelihood framework to measure the regulatory probability of genetic variants in a context-dependent manner. We show our model is superior to existing cell type-specific methods and exhibit significant GWAS signal enrichment. Using phenotypically relevant epigenomes to weight GWAS SNPs, we discover more disease-associated genes owing to regulatory changes and improve the statistical power in gene-based association test.

# Installing cepip

We have integrated our JAVA version cepip into KGGSeq package (http://grass.cgs.hku.hk/limx/kggseq/doc10/UserManual.html#ContexSpcific), now it supports whole genome scoring in an efficient manner (require mannually download dbNCFP wohle genome annotation file). Please refer to the manual for package configuration and running parameters (--regulatory-causing-predict and --cell). 
Old infomation about installation, usage and updates, please visit: http://mulinlab.org/cepip
