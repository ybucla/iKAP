# iKAP

Introduction
------------
Protein phosphorylation is one of the most important post-translational modifications (PTMs) in living cells, and is involved in the regulation of a large number of biological process, such as cell cycle, signal transduction, differentiation, proliferation, metabolism, and autophagy (Manning et al., 2002; Xie et al., 2015). With the development of liquid chromatography-tandem mass spectrometry (LC/MS-MS), proteomic technology provides the possibility to study the phosphorylation events in detail (Huttlin et al., 2010;Mann et al., 2002). Currently, Several phosphoproteomic studies were performed to characterize the landscape of phosphorylation in various system (Yeh et al., 2011;Bennetzen et al., 2012). However, the identification of key upstream regulators, namely protein kinases, from the phosphoproteomic data is still difficult. In this study, we employed stable isotope labeling by amino acids in cell culture (SILAC) technique and systematically quantified the phosphoproteome in a mouse neuroblastoma cell line N2a with treatments of Cory and Cory B. We developed a novel algorithm of in silico Kinome Activity Profiling (iKAP) to computationally infer the kinase activities from the phosphoproteomic data.

Overview
------------
The iKAP 1.0 workflow consists of: 
(1) (Optional) Pre process with SILAC based Phosphoproteomics data identified with MaxQuant and generate ELM file 
(2) Kinase prediction using iGPS software. 
(3) KA analysis.