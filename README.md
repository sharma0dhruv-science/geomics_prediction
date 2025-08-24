# 🧬 Genomics Variant Effect Prediction

This project uses **machine learning** to predict whether genetic variants (mutations) are likely to be **benign** (harmless) or **pathogenic** (disease-causing).  

---

## 🌍 Why This Project Matters
Every person’s DNA carries thousands of variants, but only a small fraction may lead to disease.  
Determining which variants are harmful is a key challenge in **genomics, precision medicine, and genetic counseling**.  
This project demonstrates how computational approaches can help classify variants quickly and consistently.

---

## 🚀 How It Works

### 1. Data Input
- Uses **ClinVar**, a public database of human genetic variants.  
- Each row contains:
  - Chromosome & position (`chrom`, `pos`)  
  - DNA change (`ref`, `alt`)  
  - Protein change (`hgvs_p`)  
  - Expert label: **Pathogenic** or **Benign**

---

### 2. Feature Engineering
To make variants “machine-readable,” biological properties are converted into **numeric features**, such as:
- **Grantham score** – chemical severity of the amino acid substitution  
- **Conservation** – whether the site is preserved across evolution (important positions are more sensitive to mutations)  
- **Protein domain info** – whether the change occurs in a functional domain  
- **Charge & hydrophobicity changes** – shifts in molecular properties

---

### 3. Training the Model
Two ML algorithms are applied:
- **Logistic Regression** – interpretable baseline  
- **Random Forest** – ensemble model that captures nonlinear effects

Data is split into:
- **Training set** – variants the model learns from  
- **Test set** – new variants used to check accuracy

---

### 4. Evaluation
The models are compared using:
- **ROC-AUC & PR-AUC curves** – ability to distinguish pathogenic from benign variants  
- **Confusion matrix** – correct vs incorrect predictions  
- **Feature importance** – which biological features matter most

---

### 5. Prediction
The best model is saved and can predict the probability of a **new mutation** being pathogenic.  

Example usage:
```python
features = {
    "grantham": 120,
    "is_conserved": 1,
    "in_domain": 1,
    "hydrophobicity_delta": 0.2,
    "charge_delta": 1
}
prob = model.predict_proba([features])[:,1]
print("Probability pathogenic:", prob)

### To get the datasets
curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Run the sampling command to trim the dataset
zcat variant_summary.txt.gz | head -10000 > variant_summary_sample.txt

Use the converter to prep the dataset
python scripts/clinvar_to_csv.py --input data/variant_summary_sample.txt --output data/clinvar_clean.csv




