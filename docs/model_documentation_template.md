This is a comprehensive **Model Documentation Template** that consolidates the structure, checklists, assumption registers, and governance text we discussed.

It is designed to be "fill-in-the-blank" compliant with **PRA SS1/23** (UK) and **SR 11-7** (US).

---

# Model Documentation Template

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)
**Standard:** Replicability & Independent Validation

---

## Document Control & Metadata

*Governance ownership tracking (SS1/23 Principle 2).*

| Field | Details |
| --- | --- |
| **Model Name** | `[e.g., IFRS9 Retail Credit Risk Scorecard]` |
| **Model ID** | `[Internal Inventory ID]` |
| **Model Tier / Materiality** | `[e.g., Tier 1 (Critical)]` |
| **Model Owner** | `[Accountable Executive]` |
| **Model Developer(s)** | `[Technical Team Names]` |
| **Version** | `[e.g., v2.1]` |
| **Last Approval Date** | `[Date]` |
| **Review Frequency** | `[e.g., Annual]` |

---

## 1. Executive Summary

*High-level synopsis for Senior Management.*

* **Purpose:** What business decision does this model support?
* **Methodology:** One-sentence summary (e.g., "Logistic regression using macroeconomic variables").
* **Performance:** Key metrics (e.g., "Gini: 65%, Stable").
* **Primary Limitation:** The single most critical risk or weakness.

---

## 2. Model Scope and Usage

*Defined scope to prevent "Model Misuse" (SS1/23 Principle 1).*

### 2.1 Business Purpose

* **Problem Statement:** The specific business problem being solved.
* **Intended Use:** How outputs are applied (e.g., "Input to Capital Calculation").
* **User Group:** Who consumes the outputs?

### 2.2 Scope and Coverage

* **Target Population:** (e.g., "UK Mortgages originated > 2015").
* **Exclusions:** Explicitly list what is out of scope.

---

## 3. Mathematical Theory and Methodology

*Justification of the approach against alternatives (SR 11-7).*

### 3.1 Theoretical Framework

* **Methodology:**
* *Example:* "Cox Proportional Hazards model defined as:"




* **Justification:** Why this method? (e.g., "Interpretability required by regulation").

### 3.2 Alternatives Considered

* **Alternative 1:** (e.g., Random Forest).
* **Reason for Rejection:** (e.g., "Lack of transparency/explainability").

### 3.3 Key Assumptions Register

*Evaluation of assumptions and their materiality.*

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | *Linearity* | Tested via Box-Tidwell. | **High:** Underestimation of risk at extremes. | Bivariate plots. | **Amber** |
| A.2 | *Economic Stability* | Historical data covers full cycle. | **Critical:** Parameters invalid in stress. | Quarterly stress test. | **Red** |

---

## 4. Data

*Data lineage and quality (SS1/23 Principle 3).*

### 4.1 Data Sources

* **Source Systems:** Database names/Schema.
* **Extraction Logic:** Reference to SQL/Python scripts.

### 4.2 Processing & Cleaning

* **Imputation:** Handling of missing values (e.g., Mean imputation vs. Exclusion).
* **Outlier Treatment:** (e.g., Winsorization at 99th percentile).

### 4.3 Data Quality & Proxies

* **Quality Checks:** Summary of completeness/accuracy checks.
* **Proxies:** Justification for any proxy data used (e.g., "Euribor used for Funding Cost").

---

## 5. Model Development

*The "Recipe" for replication.*

### 5.1 Variable Selection Checklist

*Evidence of unbiased selection.*

* **Univariate Filter:**
* [ ] Removed variables with >20% missing data.
* [ ] Removed variables with near-zero variance.
* [ ] Verified logical sign/direction against target.


* **Statistical Filter:**
* [ ] Information Value (IV) threshold applied.
* [ ] Multicollinearity check (VIF < 5).
* [ ] Correlation check (Pearson < 0.7).


* **Regulatory Filter:**
* [ ] **Protected Characteristics:** Verified no usage of race, gender, religion (UK Equality Act 2010).



### 5.2 Parameter Estimation

* **Training Data:** Definition of the training set (Time period, Split %).
* **Optimization:** Method used (e.g., Maximum Likelihood).
* **Final Equation:** The explicit model specification.

---

## 6. Testing and Performance

*Evidence of fitness for purpose.*

### 6.1 Statistical Performance

* **In-Sample:** Fit metrics (AIC, R-Squared).
* **Out-of-Sample:** Validation metrics (RMSE, AUC, Confusion Matrix).
* **Stability:** PSI (Population Stability Index) results.

### 6.2 Stress Testing & Sensitivity

* **Sensitivity:** Impact of +/- 10% change in inputs.
* **Stress Scenarios:** Performance under severe economic conditions (linked to "Red" assumptions).

### 6.3 Benchmarking

* Comparison against a "Challenger Model" or naive baseline.

---

## 7. Limitations and Post-Model Adjustments (PMAs)

*Management of uncertainty.*

### 7.1 Known Limitations

* Scenarios where the model performs poorly.
* **Mitigation:** Operational controls or buffers applied.

### 7.2 Post-Model Adjustments (Overlays)

* **Rationale:** Why is the adjustment needed?
* **Calculation:** How is the value derived?
* **Approval:** Governance body sign-off.

---

## 8. Ongoing Monitoring and Governance

*SS1/23 "Three Lines of Defense" structure.*

### 8.1 Governance Framework

> This model operates under the Bank's Model Risk Management (MRM) Framework.
> * **1st Line (Developer):** Daily performance monitoring.
> * **2nd Line (Validation):** Independent challenge and periodic review.
> * **3rd Line (Audit):** Assurance of framework compliance.
> 
> 

### 8.2 Performance Triggers

*Action plan for performance deterioration.*

| Metric | Green (No Action) | Amber (Watch List) | Red (Remediation) |
| --- | --- | --- | --- |
| **Stability (PSI)** |  |  |  |
| **Accuracy (Gini)** | Variation  | Variation  | Variation  |

---

## 9. Appendices

* **A. Data Dictionary** (Full variable definitions)
* **B. Full Regression Output** (Coefficients, P-values, SE)
* **C. Code Snippets** (Core Logic)
* **D. References** (Regulations/Papers)

---

### Final Check for the Author

Before submitting, ensure:

1. **Replicability:** Can a colleague rebuild the model using *only* this document?
2. **Version Control:** Does the doc version match the code version?
3. **Independence:** Is the tone objective and critical, not "selling" the model?

**Would you like me to create the corresponding "Validation Testing Plan" template that an independent validator would use to review this document?**