
# 02_statistics

This code depends on the `sympy`, `numpy`, `matplotlib`, and `scipy` libraries. Make sure they are installed on your system. If not, please run the following:

```bash
pip install -r requirements.txt
```
Alternatively, you can install the required packages individually:

```bash
pip install numpy matplotlib scipy sympy
```
---

## Overview

This program simulates the empirical distribution of the **Pearson correlation coefficient** \(r\) under the null hypothesis that two variables are independent and uniformly distributed.

Specifically, for given sample sizes \(n\) (e.g., 4, 5, 6), it repeatedly generates pairs of independent uniform random samples and computes the Pearson correlation coefficient for each pair. The resulting distribution of \(r\) values approximates the sampling distribution of the correlation coefficient when there is no true correlation.

The code then plots histograms of these empirical distributions for different sample sizes, allowing visualization of how the distribution of Pearson's \(r\) changes with sample size under independence.

---

## How to use

- Run the script with `python3 prob2.py`.
- It will generate and save a plot named `pearson_distribution.png` in the working directory.
- The plot shows histograms of simulated Pearson \(r\) for sample sizes 4, 5, and 6 with 10,000 replications each.

---
