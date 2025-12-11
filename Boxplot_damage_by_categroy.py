import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('genetic_load_with_breed.csv', sep=None, engine='python')

df.columns = df.columns.str.strip()

print(df.columns)  # just to check once

plt.figure(figsize=(8, 6))

df.boxplot(column='DamagingAltAlleleCount', by='Category')

plt.xlabel('Category')
plt.ylabel('Number of damaging alternate alleles')
plt.title('Genetic load by categories')
plt.suptitle('')

plt.tight_layout()
plt.show()  