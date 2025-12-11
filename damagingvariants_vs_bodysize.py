import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('genetic_load_with_breed.csv', sep=None, engine='python')

df.columns = df.columns.str.strip()
dogs_only = df[df['Category'] == 'Breed_Dogs']


print(df.columns) 

plt.figure(figsize=(8, 6))

plt.scatter(dogs_only['weight_value'], dogs_only['DamagingAltAlleleCount'], alpha=0.5)

plt.xlabel('Body weight (kg)')
plt.ylabel('Number of damaging alternate alleles')
plt.title('Genetic load vs body size')

plt.tight_layout()
plt.show()