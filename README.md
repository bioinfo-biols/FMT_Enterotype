## The interplay of gut microbiota between donors and recipients determines the efficacy of fecal microbiota transplantation
### Our main findings include: 
* the congruence of gut microbiota dysbiosis shared by CDI and IBD can be used to characterize and differentiate patients, who are classified into two common enterotypes, RCPT/E (dark red) and RCPT/B (orange).
* recipients with different enterotypes exhibited strong donor preference in FMT, with DONOR/P being more suitable for RCPT/B, particularly in IBD.
* we proposed a novel measure, the ratio of colonizers to residents after FMT (C2R), to quantify the engraftment of donor-derived bacteria in recipient, which is significantly linked to FMT outcomes and scalable for metagenomic datasets.
* we constructed an enterotype-based donor selection (EDS) model for donor-recipient matching with AUROC > 0.8.

<p align="center">
  <img src=https://user-images.githubusercontent.com/34981680/157368075-a2baf268-f4f2-49de-9090-cac877f9b425.png>
</p>
<p align="center">
  Our donor selection model based on enterotype for FMT (EDS model) performance in validations.  
</p>

### Directory structure:  
&nbsp; &nbsp; &nbsp; &nbsp; **analysis_code_Rmd**: pure code of our analyses, which is suitable to read online.  
&nbsp; &nbsp; &nbsp; &nbsp; **analysis_figures_Rnotebook**: PDF and HTML versions of our analysis to make it easier to understand our analyses with the output and figures. It is best to download the raw file, and read it with a PDF reader or browser.    
&nbsp; &nbsp; &nbsp; &nbsp; **C2R**: code for quantify the ratio of colonizers to residents after FMT (C2R).  
    
### Engraftment ratio usage:
    source('C2R.R')
    result <- C2R(FMT_abundance, FMT_config, group = 'Outcomes', blocked = 0)
Input: microbial abundance matrix and FMT config file, formatted as test files  
Output: the ratio of colonizers to residents after FMT (C2R) and its p value with group (FMT outcomes)  
  
We provided test Rmd script and test files in the folder "C2R":  
&nbsp; &nbsp; &nbsp; &nbsp; test script: &nbsp; &nbsp; test_C2R.Rmd  
&nbsp; &nbsp; &nbsp; &nbsp; test files: &nbsp; FMT_abundance.tsv and FMT_abundance.tsv  

### Citation:
Ruiqiao He, Pan Li, Jinfeng Wang, Bota Cui, Faming Zhang & Fangqing Zhao (2022) The interplay of gut microbiota between donors and recipients determines the efficacy of fecal microbiota transplantation, Gut Microbes, 14:1, DOI: 10.1080/19490976.2022.2100197

