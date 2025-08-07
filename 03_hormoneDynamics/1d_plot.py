import pandas as pd
import matplotlib.pyplot as plt

# Read CSV file
data = pd.read_csv("result.csv")

# Plot
plt.plot(data["x"], data["phi"], label="phi", linestyle="dotted")
plt.plot(data["x"], data["phi_exact"], label="phi_exact", linestyle='dotted')
#plt.plot(data["x"], data["phi"]-data["phi_exact"])
plt.xlabel("x")
plt.ylabel("Phi(x)")
plt.title("After time=t_final the plot of Phi(x) vs x")
plt.grid(True)
plt.legend()
plt.show()
