# iKAP

Introduction
------------
Protein phosphorylation is one of the most important post-translational modifications (PTMs) in living cells, and is involved in the regulation of a large number of biological process, such as cell cycle, signal transduction, differentiation, proliferation, metabolism, and autophagy ([Manning et al., 2002](https://www.ncbi.nlm.nih.gov/pubmed/12471243); [Xie et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/25484070)). With the development of liquid chromatography-tandem mass spectrometry (LC/MS-MS), proteomic technology provides the possibility to study the phosphorylation events in detail ([Huttlin et al., 2010](https://www.ncbi.nlm.nih.gov/pubmed/21183079);[Mann et al., 2002](https://www.ncbi.nlm.nih.gov/pubmed/12007495)). Currently, Several phosphoproteomic studies were performed to characterize the landscape of phosphorylation in various system ([Yeh et al., 2011](https://www.ncbi.nlm.nih.gov/pubmed/21460632);[Bennetzen et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/22517431)). However, the identification of key upstream regulators, namely protein kinases, from the phosphoproteomic data is still difficult. In this study, we employed stable isotope labeling by amino acids in cell culture (SILAC) technique and systematically quantified the phosphoproteome in a mouse neuroblastoma cell line N2a with treatments of Cory and Cory B. We developed a novel algorithm of in silico Kinome Activity Profiling (iKAP) to computationally infer the kinase activities from the phosphoproteomic data.

Overview
------------
The iKAP 1.0 workflow consists of: <br />
(1) (Optional) Pre process with SILAC based Phosphoproteomics data identified with MaxQuant and generate ELM file <br />
(2) Kinase prediction using iGPS software. <br />
(3) KA analysis.<br />

Pre-requisites and Install
------------
    1. python 2.7.x.x
	2. python package 'scipy' and 'numpy' installed
	3. iGPS installed

Usage
------------
(1). (Optional) Pre process with SILAC based Phosphoproteomics data identified with MaxQuant and generate ELM file <br />
```
python getpep.py -h
Usage: python getpep.py [options] -i MSresultfile -s seqfile -u proteincolumn -p pepcolumn -r ratiocolumn -o output

Options:
  -h, --help            show this help message and exit
  -i MSFILE				quantification file using MaxQuant or others [Default none]
  -s SEQFILE			fasta file used for the database search [Default none]
  -o OUTDIR				output directory [Default output]
  -u PROTEIN            protein column, 0-based, like "IPI00021812.2" [Default 0]
  -p PEP                pep column, 0-based, like "_HRS(ph)NS(ph)FSDER_" [Default 1]
  -r RATIO              ratio column, 0-based, like "0.38957" [Default 4]
```

(2). Kinase prediction using iGPS software. <br />
```
1. Open iGPS 1.0 and COPY the peptides from the previous result "PhosphoPep.txt" into iGPS software. 
2. Select All kinases. 
3. Set species to "M.musculus" (or species in research) and thresehold as "Low". 
4. Submit and wait for the calculation, which may take about 1 min according the amount of peptides input. 
5. The prediction results will be displayed in the result panel. 
6. Save the results to file. Here, the results were saved into "output/kinase_prediction.iGPS".
```
![](https://github.com/ybucla/iKAP/blob/master/igpsflow.png)

(3). KA analysis.<br />
```
python KAnalysis.py -h
Usage: python KAnalysis.py -i ratioElmfile -g iGPSResultfile -o output

Options:
  -h, --help   show this help message and exit
  -i ELMFILE   quantification ratio elm file [Default none]
  -g IGPSFILE  result file from iGPS [Default none]
  -o OUTDIR    output directory [Default output]
```

