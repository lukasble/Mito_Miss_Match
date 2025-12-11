import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('genetic_load_with_breed.csv', sep=None, engine='python')

print(df.columns)


plt.figure(figsize=(10,6))
plt.hist(df['DamagingAltAlleleCount'], bins=40, color='steelblue', edgecolor='black')
plt.ylabel('Number of Damaging Alternate Alleles')
plt.xlabel('Number of Individuals')
plt.title('Distribution of Damaging Allele Counts Across All Individiuals')

plt.tight_layout()
plt.show()