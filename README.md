# ECAT
Evolutionary Cluster-based Association Test (ECAT) is an adaptation of linear mixed modles (LMMs) to incorporate homoplasy information for enhancing the identiﬁcation of genes associated with drug resistance. ECAT involves three phases: homoplasy inference, detection of clustered (i.e. hyper-variable) regions, and association tests using LMMs. Here are the prerequisites and scripts for running the ECAT pipeline.


- [Required packages](#required-packages)
- [Input data format](#input-data-format)
- [Running ECAT](#running-ecat)
- [Citing ECAT](#citing-ecat)
- [License](#license)

# Required packages

- Phylogenetic tree reconstruction: 
[PAUP](https://paup.phylosolutions.com)
- Association test: 
[GEMMA](https://github.com/genetics-statistics/GEMMA) 


# Input data format

After calling variants from the assembled genomes, we generate a multiple sequence alignment (MSA) of all samples by aligning them to a reference genome. Here, we use a collection of 30 clinical isolates of *Mycobacterium tuberculosis* and the reference genome H37Rv (GenBank: NC_000962.2) as an example, sample_30h.align. The format of the MSA is descibed in details as follows. The first 31 lines of the MSA consists of the orders of strains where the first one is the reference. Other lines represent the aigned sequences for sites across the genome. The sites exhibiting SNPs are labeled with asterisks and annotated. The format of information for each site are coordinate, gene name, alias of the gene name if any, a dot symbol, an aligned sequence of nucleotides, an asterisk if it is a polymorphic site , an annotation if the SNP occurs within a gene.
    
```
# H37Rv
# Peru6143h
# Peru6145h
# Peru6146h
# Peru6148h
# Peru6149h
# Peru6150h
# Peru6151h
# Peru6152h
# Peru6153h
# Peru6154h
# Peru6158h
# Peru6159h
# Peru6161h
# Peru6162h
# Peru6163h
# Peru6164h
# Peru6174h
# Peru6176h
# Peru6177h
# Peru6179h
# Peru6181h
# Peru6182h
# Peru6189h
# Peru6194h
# Peru6229h
# Peru6230h
# Peru6231h
# Peru6232h
# Peru6237h
# Peru6238h
1       Rv0001     dnaA       . TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  
2       Rv0001     dnaA       . TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  
3       Rv0001     dnaA       . GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
...
71      Rv0001     dnaA       . CCCCCCCCCTCCCCCCCCCCTCCCCCCCCCC * T:[P24L]
72      Rv0001     dnaA       . TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  
73      Rv0001     dnaA       . AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  
74      Rv0001     dnaA       . AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  
75      Rv0001     dnaA       . GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG  
76      Rv0001     dnaA       . GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG  
77      Rv0001     dnaA       . TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  
78      Rv0001     dnaA       . TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  
79      Rv0001     dnaA       . GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG  
80      Rv0001     dnaA       . AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  
81      Rv0001     dnaA       . CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
82      Rv0001     dnaA       . GGGGGGGGGGGGGGGGGGGGGCGGGGGGGGG * C:[D28H]
...
```   
 

## Data preprocessing

1. We obtain all SNPs by running

```
python get_snp.py sample_30h.align > sample_30h.snp
```

   , where we exclude ambiguous sites, repetive regions, and PE/PGRS genes. 


2. We use both synonymous and nonsynonymous SNPs to generate a phylogenetic tree by maximum parsomony using PAUP. The tree is named as sample_30h.tre.


3. We filter out synonymous SNPs by running

```
python filter_synon_snp.py sample_30h.snp > sample_30h_nonsyn.snp
```

   , where we keep all SNPs within the intergenic regions and *rrs* (16S rRNA).



# Running ECAT

  The input files are an alignment of all nonsynonymous SNPs, a phylogenetic tree, and a list of phenotypic traits of interests.


## Phase I--Homoplasy Inference

- We calculate the homoplasy count for each site by running

```
python count_HI.py sample_30h.tre sample_30h_nonsyn.snp > sample_30h_nonsyn_HI.snp
```

   , where the last column of each site is the homoplasy index (HI), which is the number of excess changes plus 1. Thus, a site with 1 HI means that it does not need extra changes after being mapped onto the tree, which is homoplasy-free. In contrast, a site with HI over 1 represents it is homoplasic.


## Phase II--Clustered Region Identification

- We identify the clustered regions based on the homoplasy indexes by running

```
python chuck_genome_cluster_locmu_HI.py sample_30h_nonsyn_HI.snp 7500 20 > sample_30h_nonsyn_HI_region_15kwin20
```

   , where the second argument, 7500, is the half size of a sliding window for calculating the local mutation rate 
and the third argument, 20, is a given span of SNPs as a region that we group adjacent sites up to the span.

  Here, we will obtain a list of non-overlapping clustered regions sorted by the adjusted *p*-values with 5% FDR cutoff.  


## Phase III--Association Test 

- Phenotypes

  The phenotypes of all strains are converted to 0/1 (S/R) and recorded in the same order as in the genotype file. The file is named as s30h_XXX_bin.txt, where XXX is the abbreviation of the drug (e.g. INH, RIF, EMB, etc.).

- Genotypes

  We apply a burden test to group SNPs within each region and convert them to the format required by GEMMA.

```
python convert_gemma_format_region.py sample_30h_nonsyn.snp sample_30h_nonsyn_HI_region_15kwin20 > sample_30h_nonsyn_HI_region_gemma_15k
```

- Genetic relatedness matrix in a LMM

  To calculate a genetic relatedness matrix (kinship), we firstly filter out known drug-resistant sites/genes by running 

```
python filter_ds.py sample_30h_nonsyn.snp > sample_30h_nonsyn_fil.snp
```

   , and then convert it to the format required by GEMMA

```
python convert_gemma_format.py sample_30h_nonsyn_fil.snp > sample_30h_fil_site_gemma
```

   and finally compute the centered kinship by running

```
./gemma.linux -g sample_30h_fil_site_gemma -p s30h_INH_bin.txt -maf 0 -gk 1  -o sample_30h_grm 
```

- Association Test using GEMMA

```
./gemma.linux -g sample_30h_nonsyn_HI_region_gemma_15k  -p s30h_INH_bin.txt -maf 0 -k $PATH/output/sample_30h_grm.cXX.txt -lmm -o sample_30h_region_HI_15k_INH
```

- Multiple test correction 

  For multiple test coorection, we adjust the Wald P values by 5% flase discovery rate (FDR), annotate the regions, and sort the results

```
python assign_pwald_region_ant_pos.py PATH/output/sample_30h_region_HI_15k_INH.assoc.txt s30h_INH_bin.txt sample_30h_nonsyn_HI_region_gemma_15k H37Rv.prot_table sample_30h_nonsyn_HI_region_15kwin20 | sort -gk 21 > sample_30h_region_homo_15k_INH_lmm_sorted 
```

 

# Citing ECAT

# License

Authors: Yi-Pin Lai and Thomas R. Ioerger.

ECAT is a free software: you can redistribute it and/or modify it under the terms of the [GNU General Public License](http://www.gnu.org/licenses/) as published by the Free Software Foundation.

