# 2025_Spatial_Profiling_in_Medullary_Extramedullary_AML

<h3 align="center">
Spatial Transcriptomic Profiling of Medullary and Extramedullary Acute Myeloid Leukemia (AML)
</h3>

<p align="center">
<strong>Dasdemir et al., 2025, iScience</strong>  
<br/>
<em>"Integrative Spatial Multi-Omics Reveal Niche-Specific Inflammatory Signaling and Differentiation Hierarchies in Acute Myeloid Leukemia"</em>
</p>

---

## 🧾 Overview

This repository contains the code used for the analysis presented in the manuscript:

> **Integrative Spatial Multi-Omics Reveal Niche-Specific Inflammatory Signaling and Differentiation Hierarchies in Acute Myeloid Leukemia**  
> *Enes Dasdemir, Ivo Veletic, Christopher P Ly, Andres E Quesada, Christopher D Pacheco, Fatima Z Jelloul, Pamella Borges, Sreyashi Basu, Sonali Jindal, Zhiqiang Wang, Alexander Lazar, Khalida M Wani, Dinler A Antunes, Patrick K Reville, Preethi H Gunaratne, Robert J Tower, Padmanee Sharma, Hussein A Abbas, 2025, iScience*

The project integrates:

- **10x Genomics Visium** spatial transcriptomics  
- **GeoMx** spatial profiling  
- Multiplex immunofluorescence / histology context  

to characterize:

- Spatial organization of leukemic cells in **medullary** and **extramedullary** tissues  
- Inflammatory microenvironmental niches  
- Trans-differentiation states and niche-specific AML programs  
- Cell–cell communication and inflammatory signaling axes

---

## 📁 Repository Structure

```text
.
├── R/Initial_QC.R
├── R/Clustering_Annotation.R
├── R/Downstream_Analysi
├─- R/MAD_Filtering_Process.R
├── Opal_Intensities
└── README

