# DM-CART: Dirichlet-Multinomial Classification and Regression Trees

## Master's Thesis
**Title:** Investigating Determinants of Birth Weight Using Bayesian Tree-Based Nonparametric Modeling  
**Author:** Adam Kurth  
**Department/Degree:** Department of Mathematical & Statistical Sciences, Arizona State University — Masters of Statistics in Statistics
**Advisor/Committee Chair:** [Dr. P. Richard Hahn](https://math.la.asu.edu/~prhahn/) 

## Overview

This repository contains the code, data, and results for my Master's thesis on birth weight modeling using a novel Bayesian nonparametric framework called *DM-CART* (Dirichlet-Multinomial Classification and Regression Trees). The project extends traditional CART methodology in a new direction by applying Dirichlet-Multinomial distributions to model the full distribution of birth weight outcomes. The goal is to provide a robust and interpretable framework for predicting low birth weight (LBW) risk based on maternal-infant characteristics. The DM-CART algorithm is implemented as an extension of the `rpart` package in R, allowing for tree-based modeling of multinomial outcomes. The model is fit using the 2021 U.S. Natality dataset, which contains over 3 million births, and the previous year's data (2020) is used to inform the Dirichlet priors. The analysis focuses on identifying high-risk maternal profiles and their respective predictors for LBW, with the aim of providing actionable insights for clinicians and public health officials.


## Abstract

Low birth weight (LBW) remains a critical public-health indicator, linked strongly with higher neonatal mortality, developmental delays, and lifelong chronic diseases. Using the 2021 U.S. Natality dataset (> 3 million births), this thesis develops a Bayesian, tree-based, nonparametric framework that models the full birth weight distribution and quantifies LBW risk.

The raw dataset is condensed into 128 mutually exclusive classes defined by seven dichotomous maternal-infant predictors and 11 birth weight categories, comprised of 10\% LBW quantile categories plus one aggregated normal weight category for added LBW granularity. Classification and Regression Trees (CART) are grown using the marginal Dirichlet-Multinomial likelihood as the splitting criterion. This criterion is equipped to handle sparse observations, with the Dirichlet hyperparameters informed by previous quantiles from the 2020 dataset to avoid "double dipping".

Employing a two-tier parametric bootstrap resampling technique, a 10,000 tree ensemble is grown yielding highly stable prediction estimates. Maternal race, smoking status, and marital status consistently drive the initial LBW risk stratification, identifying Black, smoking, unmarried mothers among the highest-risk subgroups. When the analysis is restricted to LBW births only, infant sex and maternal age supersede smoking and marital status as key discriminators, revealing finer biological gradients of risk. Ensemble predictions are well calibrated, and 95\% bootstrap confidence intervals achieve nominal coverage.

The resulting framework combines the interpretability of decision trees with Bayesian uncertainty quantification, delivering actionable, clinically relevant insights for targeting maternal-health interventions among the most vulnerable subpopulations.

## Repository Structure

### Core Implementation Files

- **R Scripts**
  - `dm-cart.R`: Implementation of the DM-CART algorithm with extensions to the `rpart` package, including:
    - Dirichlet-Multinomial Marginal Likelihood splitting criterion
    - Custom initialization, evaluation, and prediction functions
    - Tree pruning and visualization 
    - Bootstrap resampling for uncertainty quantification
    - Variable importance analysis
 
  - `data_rebin.R`: Data preprocessing, creation of maternal-infant profiles and multinomial counts data ($\mathbf{X}$ and $\mathbf{Y}$, respectively)

  - `quantile.R`: Creating quantile-based birth weight categories and cutpoints for consistent binning

### Data Organization

- **Data Directory**
  - `birthweight_data/`: Contains the original U.S. Natality datasets and processed data files
    - `natality2020us-original.csv`: Raw 2020 data to obtain Dirichlet prior quantiles 
    - `natality2021us-original.csv`: Raw 2021 data used in analysis
    - `rebin/`: Processed data with 2.5kg threshold included 
      - `rebin_with_2.5kg/`: Full model data, processed *including* 2.5kg +
      - `rebin_without_2.5kg/`: LBW-only model data, processed *excluding* 2.5kg +
    - $\mathbf{X}$: fixed predictor matrix for maternal-infant profiles (128 rows × 7 columns)
    - $\mathbf{Y}$ : consolidated counts matrix across 128 profile classes and 10 or 11 birth weight categories (128 rows × 10/11 columns).
## Methodology

The DM-CART algorithm is designed to model the distribution of birth weights using a tree-based approach.

### 1. Data Preprocessing and Discretization
  
  - Birth weights $\leq 2.5$ kg (LBW) are allocated across 10 quantile-based categories, called LBW-region
  - Full model includes all observations above 2.5kg threshold (i.e. normal birth weight (NBW) (>2500g)) while the LBW-only model does not
  
- **Data Preparation**:
  - Seven maternal-infant characteristics converted to binary predictors
  - Creates 128 ($2^7$) distinct maternal-infant profiles
  - Counts for each profile and birth weight category form a multinomial contingency table, called counts data

### 2. DM-CART Methodology

- **Theoretical Foundation**:
  - Extension of CART methodology to handle multinomial counts data for birth weight modeling
  - Utilizes the Dirichlet-Multinomial distribution for counts modeling
  - Dirichlet prior provides natural zero-count handling and smoothing for sparse data
  - Uses adjusted marginal likelihood as the splitting criterion and objective function

- **Model Components**:
  - **Custom `rpart` Method**: Extends R's `rpart` package with specialized functions
  - **Splitting Criterion**: Uses the log-likelihood difference between splitting versus not splitting
  - **Evaluation Method**: Computes Dirichlet-Multinomial log-likelihood at each node
  - **Prediction Function**: Generates probability distributions with Dirichlet smoothing

- **Implementation Details**:
  ```r
  # Implementation structure in custom-rpart.R:
  myinit <- function(y, offset, parms, wt) {...}  # Initialize DM tree
  myeval <- function(y, wt, parms) {...}          # Evaluate nodes with DM likelihood
  mysplit <- function(y, wt, x, parms, continuous) {...}  # DM likelihood ratio splitting
  mypred <- function(fit, newdata, type) {...}    # Prediction with Dirichlet smoothing
  dm.method <- list(init=myinit, eval=myeval, split=mysplit, pred=mypred, method="dm")
  ```

## Installation and Usage

### Prerequisites

- R (version 4.0 or higher)
- Required R packages:
  - `rpart`: Foundation for tree-based models (with custom extensions)
  - `ggplot2`, `reshape2`, `dplyr`: For data visualization and manipulation
  - `xtable`: For LaTeX table generation

### Data Requirements

1. **U.S. Natality datasets** for 2020 and 2021:
   - Available from the National Center for Health Statistics (NCHS)
   - [NCHS Natality Data](https://www.cdc.gov/nchs/data_access/Vitalstatsonline.htm)
   - Download the datasets and place them in the `birthweight_data/` directory
   - Note: The datasets are large (>3 GB each) and not included in this repository due to size constraints
   - Requires variables:
     - `dbwt`: Birth weight in grams (response variable)
     - `sex`, `dmar`, `mrace15`, `mager`, `meduc`, `precare5`, `cig_0`: Predictors

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/dm-cart.git
   cd dm-cart
   ```

2. Install required R packages:
   ```r
   install.packages(c("rpart", "ggplot2", "xtable", "rnn", "reshape2", "dplyr"))
   ```

3. Download and place the U.S. Natality datasets in the `birthweight_data/` directory
   - Files should be named: `natality2020us-original.csv` and `natality2021us-original.csv`

4. Run the data preprocessing script:
   ```r
    source("data_rebin.R")
    ```
    - This will create the processed data files in the `birthweight_data/rebin/` directory 

5. Run the DM-CART algorithm:
   ```r
   source("dm-cart.R")
   ```
   - This will fit the DM-CART model to the processed data and generate bootstrap ensemble
  
### Key Customization Options

- **Tree Complexity**: Modify `rpart.control()` parameters:
  ```r
  dm.control <- rpart.control(minsplit=2, cp=0, maxdepth=8, xval=0)
  ```

- **Bootstrap Size**: Adjust the number of bootstrap iterations:
  ```r
  B <- 10000  # Default is 10,000 iterations
  ```

## Data Description

### Predictor Variables

The model uses seven binary predictors that define 128 (2^7) distinct maternal-infant profiles:

| Variable | Source | Description | Encoding |
|----------|--------|-------------|----------|
| **BOY** | `sex` | Infant sex | 1 = Male, 0 = Female |
| **MARRIED** | `dmar` | Maternal marital status | 1 = Married, 0 = Not married |
| **BLACK** | `mrace15` | Maternal race | 1 = Black only (mrace15=2), 0 = Other |
| **OVER33** | `mager` | Maternal age | 1 = Over 33 years, 0 = 33 years or younger |
| **HIGH SCHOOL** | `meduc` | Maternal education | 1 = High school graduate (meduc=3), 0 = Other education level |
| **FULL PRENATAL** | `precare5` | Prenatal care received | 1 = Full prenatal care (precare5=1), 0 = Less than full care |
| **SMOKER** | `cig_0` | Smoking during pregnancy | 1 = Any smoking (cig_0>0), 0 = No smoking |

### Response Variable
All birth weight (in g) observations are allocated into quantile-based categories (Type 1/Type 2) based on the following quantile ranges during the preprocessing procedure. The quantile ranges are defined as follows:

| **Quantile** | **Range (g)** (Type 1: LBW + Normal) | **Prior (%)** | **Range (g)** (Type 2: LBW only) | **Prior (%)** |
|--------------|--------------------------------------|---------------|----------------------------------|----------------|
| Q1           | 227–1170                             | 0.84          | 227–1170                         | 10             |
| Q2           | 1170–1644                            | 0.84          | 1170–1644                        | 10             |
| Q3           | 1644–1899                            | 0.83          | 1644–1899                        | 10             |
| Q4           | 1899–2069                            | 0.83          | 1899–2069                        | 10             |
| Q5           | 2069–2183                            | 0.87          | 2069–2183                        | 10             |
| Q6           | 2183–2270                            | 0.83          | 2183–2270                        | 10             |
| Q7           | 2270–2350                            | 0.86          | 2270–2350                        | 10             |
| Q8           | 2350–2410                            | 0.93          | 2350–2410                        | 10             |
| Q9           | 2410–2460                            | 0.71          | 2410–2460                        | 10             |
| Q10          | 2460–2500                            | 0.80          | 2460–2500                        | 10             |
| Normal       | >2500                                | 91.67         |                                  |                |


## Acknowledgements

- Thesis advisor: [Dr. P. Richard Hahn](https://math.la.asu.edu/~prhahn/), Department of Mathematical & Statistical Sciences, Arizona State University
- Committee members: 
  - Dr. [Shuang Zhou](https://sites.google.com/view/shuangzhousomss), Department of Mathematical & Statistical Sciences, Arizona State University
  - Dr. [Shiwei Lan](https://math.la.asu.edu/~slan/), Department of Mathematical & Statistical Sciences, Arizona State University
- Special thanks to the [Department of Mathematical & Statistical Sciences](https://math.asu.edu) at Arizona State University for their support

## Citation

If you use this code or methodology in your research, please cite:

```
Kurth, A. (2025). Investigating Determinants of Birth Weight Using Bayesian Tree-Based Nonparametric Modeling.
Master's Thesis, Arizona State University.
```

## Contact

For questions or further information, please contact adam.kurth@asu.edu or visit my [GitHub profile](https://github.com/adamkurth)
