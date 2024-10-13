#########################################################
# Author: Rafat Eissa 
# Date: 13th Oct 2024
# Project: Bread Wheat Assesment 
#########################################################

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.svm import SVR
from sklearn.linear_model import Ridge, Lasso
from sklearn.tree import DecisionTreeRegressor
import matplotlib.pyplot as plt


##################################################################################################################
################ 1. Read the data 
##################################################################################################################

geno_data = pd.read_csv('50k_genotype_Bread.csv')  
pheno_data = pd.read_csv('phenotype_Bread_for_Tassel.csv') 

##################################################################################################################
################ 2. Preprocess the data
##################################################################################################################

# Handling missing values (mean imputation for genotypic data)
geno_data.fillna(geno_data.mean(), inplace=True)

# Minor Allele Frequency (MAF) filtering (0.05 as a threshold)
def maf_filter(geno_matrix, threshold=0.05):
    allele_freq = geno_matrix.mean(axis=0) / 2
    maf = np.minimum(allele_freq, 1 - allele_freq)
    return geno_matrix.loc[:, maf >= threshold]

geno_data_filtered = maf_filter(geno_data)

# (OPTIONAL) Build kinship matrix (to be used for later analyses, like GWAS)
kinship_matrix = np.dot(geno_data_filtered, geno_data_filtered.T)

# (OPTIONAL) PCA for dimensionality reduction
scaler = StandardScaler()
geno_scaled = scaler.fit_transform(geno_data_filtered)
pca = PCA(n_components=10)  # Use 10 components for PCA
geno_pca = pca.fit_transform(geno_scaled)



##################################################################################################################
################ 3. Build models with different algorithms
##################################################################################################################
# Ridge Regression, Lasso Regression, SVR, Decision Trees
models = {
    "SVR": SVR(kernel='rbf'),
    "Ridge": Ridge(alpha=1.0),
    "Lasso": Lasso(alpha=0.1),
    "DecisionTree": DecisionTreeRegressor(max_depth=5)
}

X = geno_pca
y = pheno_data.values.ravel()


##################################################################################################################
################ 4. Cross-validation (K-fold)
##################################################################################################################
kf = KFold(n_splits=5, shuffle=True, random_state=42)
best_model = None
best_score = -np.inf
model_scores = {}

for name, model in models.items():
    cv_scores = cross_val_score(model, X, y, cv=kf, scoring='r2')
    mean_cv_score = np.mean(cv_scores)
    model_scores[name] = mean_cv_score
    print(f"Model: {name}, Mean CV R2 Score: {mean_cv_score}")
    if mean_cv_score > best_score:
        best_score = mean_cv_score
        best_model = model
##################################################################################################################
################ 5. Fit the best model to predict the highest individuals
##################################################################################################################


print(f"\nBest model is {best_model} with R2 score: {best_score}")
best_model.fit(X, y)

# Predict the genomic estimated breeding values (GEBVs)
predictions = best_model.predict(X)

# Select the top individuals based on their predicted GEBVs
top_individuals = np.argsort(predictions)[-10:]  # Top 10 individuals


plt.figure(figsize=(10, 6))
plt.scatter(geno_pca[:, 0], geno_pca[:, 1], c=predictions, cmap='viridis')
plt.colorbar(label='Predicted GEBVs')
plt.scatter(geno_pca[top_individuals, 0], geno_pca[top_individuals, 1], color='red', label='Top 10 Individuals')
plt.title('PCA of Genotypic Data with Predicted GEBVs')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.legend()
plt.show()

plt.figure(figsize=(8, 5))
plt.bar(model_scores.keys(), model_scores.values(), color='blue')
plt.title('Cross-Validation Scores for Different Models')
plt.ylabel('Mean R2 Score')
plt.show()

