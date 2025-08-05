import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def simulate_pearson_distribution(n, r_reps=10000):
    """
    For given sample size n, repeat r_reps times:
    - generate n iid 2D points (independent uniform)
    - calculate Pearson correlation coefficient
    Returns array of r coefficients.
    """
    r_values = []
    for _ in range(r_reps):
        x = np.random.rand(n)
        y = np.random.rand(n)
        r_sample, _ = pearsonr(x, y)  # sample Pearson r
        r_values.append(r_sample)
    return np.array(r_values)

# Run for various n
ns = [4, 5, 6, 7, 8, 9]
r_reps = 10000
plt.figure(figsize=(12, 8))

for n in ns:
    r_vals = simulate_pearson_distribution(n, r_reps)
    plt.hist(r_vals, bins=50, alpha=0.5, label=f'n={n}', density=True)

plt.title('Empirical distribution of Pearson correlation under null (independent data)') 
plt.xlabel('Pearson correlation coefficient r')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
plt.savefig('pearson_distribution.png')

