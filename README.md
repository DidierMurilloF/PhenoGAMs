# Smoothers in the Field: GAMs for Capturing Treatment Effects and Spatial Trends in Plant Breeding

![PhenoGAMs Cover](https://raw.githubusercontent.com/DidierMurilloF/PhenoGams/main/images/Smoothers.png)

## Overview

This repository accompanies the article **"Smoothers in the Field"**, which explores how Generalized Additive Models (GAMs) can be applied to field trial data to capture spatial patterns and estimate genotype performance. The analysis compares GAMs to classical mixed models fitted via ASReml, using real p-rep trial data.

📘 **Article preview**: [https://didiermurillof.github.io/PhenoGAMs](https://didiermurillof.github.io/PhenoGAMs)

---

## Key Features

- 🌾 Field trial data (p-rep design)
- 📈 Spatial modeling using `mgcv::bam()` and `asreml`
- 📊 Diagnostic plots, residual heatmaps, and prediction comparisons
- 🔁 Full reproducible R Markdown workflow

---

## Contents

- `index.Rmd` — Main R Markdown source
- `index.html` — Rendered HTML article
- `FARGO_pREP_2025-04-18.csv` — Sample dataset
- `gam_phenotypic.R` — Script to reproduce GAM model
- `README.md` — Project overview

