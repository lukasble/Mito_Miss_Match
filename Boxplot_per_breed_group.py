import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('genetic_load_with_breed.csv', sep=None, engine='python')
df.columns = df.columns.str.strip()

def assign_size(weight):
    if weight < 5: 
        return 'Toy'
    elif weight < 15:
        return 'Small'
    elif weight < 30:
        return 'Medium'
    elif weight < 45:
        return 'Large'
    else:
        return 'Giant'

df['SizeClass'] = df['weight_value'].apply(assign_size)


dogs_only = df[df['Category'] == 'Breed_Dogs']

plt.figure(figsize=(8, 6))


dogs_only.boxplot(column='DamagingAltAlleleCount', by='SizeClass', rot=90)

plt.xlabel('Size class')
plt.ylabel('Number of damaging alternate alleles')
plt.title('Genetic load across size classes')

plt.suptitle('')
plt.tight_layout()
plt.show()