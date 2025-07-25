# Aloe arborescens Antifungal Efficacy Study
Statistical analysis of Aloe arborescens compounds on Candida albicans morphogenesis using count data modeling and multiple hypothesis testing.

This project investigates the antifungal properties of compounds derived from *Aloe arborescens* on *Candida albicans*, a common opportunistic fungal pathogen. Using descriptive statistics and a quasi-Poisson regression model, we evaluate how different treatments affect hyphae formation — a key indicator of fungal virulence.

## 🔬 Objective
To assess whether isolated and combined Aloe-based treatments significantly reduce hyphal and morphological transitions in *C. albicans* cells, compared to a control group.

## 🧪 Dataset
The dataset includes morphological counts from 1,000 *C. albicans* cells under 7 treatment conditions (including a control group), each replicated 5 times.

### Morphological Response Types:
- No Morphology Change
- Germ Tubes Only (GT)
- Pseudohypha Only (PH)
- Buds Only
- Combinations (e.g., GT + PH)

### Treatment Groups:
1. Control (No treatment)
2. Whole Aloe Extract
3. Compound A (Isolated from Aloe)
4. Compound B (Isolated from Aloe)
5. A-SA (Sigma-Aldrich equivalent of A)
6. B-SA (Sigma-Aldrich equivalent of B)
7. A + B and A-SA + B-SA combinations

## 🛠️ Methods
- **Data Cleaning**: Renamed columns, created composite metrics (e.g., `Total_Hyphae`)
- **Visualizations**: Heatmaps, stacked bar charts, Q-Q plots, and Cook’s distance
- **Modeling**: Generalized Linear Model using Quasi-Poisson distribution
- **Hypothesis Testing**: Dunnett-style multiple comparisons with adjusted p-values

## 📈 Key Findings
- **Whole Aloe Extract** showed the greatest reduction in hyphae formation.
- **Combined Treatments (A + B or A-SA + B-SA)** also significantly reduced fungal morphogenesis.
- The model passed diagnostics including dispersion, Q-Q normality, and residual checks.

## 📂 Files Included
- R file with full analysis
- Rendered pdf report
- Morphological count dataset
- Client's presentation document to give an overview of project

## 📊 Output Snapshot
> Heatmaps, residual plots, and summary tables provided in the report.

## 📌 Conclusion
The findings support the potential of Aloe-derived treatments, especially whole extract and compound combinations, as promising antifungal agents that reduce morphological transitions in *C. albicans*.

---

Check out the full analysis in the `.html` report or explore the `Rmd` for reproducible code.
