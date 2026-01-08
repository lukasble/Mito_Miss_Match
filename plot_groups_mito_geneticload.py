import pandas as pd
import matplotlib.pyplot as plt
import os

# -------------------------
#Loading data
# -------------------------
path = "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/Results/9_deleterious_variants/mito_deleterious_counts.tsv"
df = pd.read_csv(path, sep="\t")

# -------------------------
# Defining output directory
# -------------------------
out_dir = "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/Sarina/Results/9_deleterious_variants"
os.makedirs(out_dir, exist_ok=True)

# -------------------------
# Definig the groups/categories
# -------------------------

# Coyotes to exclude
coyotes = {
    "CLATUS000001", "CLATUS000002", "CLATUS000003", "CLATUS000004"
}

# Wolves
wolves = {
"CLUPAZ000001","CLUPCN000001","CLUPCN000002","CLUPCN000003","CLUPCN000004",
"CLUPCN000005","CLUPCN000006","CLUPCN000007","CLUPCN000009","CLUPCN000010",
"CLUPEA000001","CLUPEU000002","CLUPGR000001","CLUPGR000002","CLUPGR000003",
"CLUPGR000004","CLUPGR000005","CLUPGR000006","CLUPGR000007","CLUPGR000008",
"CLUPGR000009","CLUPGR000010","CLUPGR000011","CLUPGR000012","CLUPIR000001",
"CLUPIR000002","CLUPIR000003","CLUPIR000004","CLUPIR000005","CLUPIR000006",
"CLUPKG000001","CLUPKZ000002","CLUPPT000001","CLUPPT000002","CLUPRU000001",
"CLUPRU000002","CLUPRU000003","CLUPRU000004","CLUPRU000005","CLUPRU000006",
"CLUPRU000007","CLUPRU000008","CLUPRU000009","CLUPRU000010","CLUPRU000011",
"CLUPRU000012","CLUPRU000013","CLUPRU000014","CLUPRU000019","CLUPRU000020",
"CLUPSE000001","CLUPSE000002","CLUPSE000003","CLUPSE000004","CLUPSE000005",
"CLUPSE000006","CLUPTJ000001","CLUPTJ000002","CLUPTJ000003","CLUPTJ000004",
"CLUPTJ000005","CLUPTJ000006","CLUPTJ000007"
}

# -------------------------
# Assigning a category
# -------------------------
def assign_group(sample_id):
    if sample_id in coyotes:
        return None
    if sample_id in wolves:
        return "Wolves"
    if sample_id.startswith("VILL"):
        return "Village dogs"
    return "Breed dogs"

df["group"] = df["dog_id"].apply(assign_group)
df = df.dropna(subset=["group"])

# -------------------------
# Preparing the data
# -------------------------
groups = ["Breed dogs", "Village dogs", "Wolves"]
data = [df[df["group"] == g]["deleterious_count"] for g in groups]

# -------------------------
# Plotting the data
# -------------------------
plt.figure(figsize=(6, 5))
plt.boxplot(data, labels=groups)
plt.ylabel("Number of deleterious mitochondrial variants")
plt.title("Mitochondrial genetic load across groups")
plt.tight_layout()

out_file = os.path.join(out_dir, "mito_deleterious_boxplots.png")
plt.savefig(out_file, dpi=300)

