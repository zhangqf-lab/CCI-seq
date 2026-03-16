# :sparkles: CCI-seq :sparkles: 

This repository provides the official implementation of the data analysis pipeline described in our study：

## :dart: Resolving cell-cell interaction networks and their molecular logic in complex tissues

<p align="justify" >
Small molecules can bind RNAs to regulate their fate and functions, providing promising opportunities for treating human diseases. However, current tools for predicting small molecule-RNA interactions (SRIs) require prior knowledge of RNA tertiary structures, limiting their utility in drug discovery. Here, we present SMRTnet, a deep learning method to predict SRIs based on RNA secondary structure. By integrating <b>two large language models</b>, <b>convolutional neural networks</b>, <b>graph attention networks</b>, and <b>an attention-based multimodal data fusion model</b>, SMRTnet achieves high performance across multiple experimental benchmarks, substantially outperforming existing state-of-the-art tools. 
</p>
<p align="justify" >
For wet-lab validation, we conducted a large-scale experimental assessment on SMRTnet predictions for 10 disease-associated RNA targets (<i>e.g.</i> <b>mRNA of undruggable proteins, onco-miRNAs, viral RNAs, and RNA repeat expansions</b>), identifying 40 hits of RNA-targeting small molecules with nanomolar-to-micromolar dissociation constants using microscale thermophoresis (MST). Focusing on the <i>MYC</i> internal ribosome entry site (IRES) as a target, SMRTnet-predicted small molecules showed binding scores correlated closely with observed validation rates. Notably, one predicted compound downregulated <i>MYC</i> expression, inhibited proliferation, and promoted apoptosis in three cancer cell lines. 
</p>
<p align="justify" >
Taken together, SMRTnet expands the scope of feasible RNA targets and accelerates the discovery and development of RNA-targeting therapeutics.
</p>

<p align="center"><img src="figs/workflow.png" width=100% /></p>
<p align="center" > <b>Overview of SMRTnet</b> </p>




## Resolving cell-cell interaction networks and their molecular logic in complex tissues

## Description
Jupyter notebooks for CCI-seq data analysis


Processed files are available for download at https://doi.org/10.5281/zenodo.17445811


<!-- CONTACT -->
## Contact
Kang Tian - tiankang@mail.tsinghua.edu.cn
<p align="right">(<a href="#readme-top">back to top</a>)</p>
