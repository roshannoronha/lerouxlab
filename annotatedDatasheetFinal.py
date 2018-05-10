#Written by Roshan Noronha
#April 4, 2018

import pandas as pd
import numpy as np

#read in Leroux/Ritter spreadsheet
lr_csv = pd.read_csv('lerouxRitter.csv')
#read in human ciliary genes spreadsheet
hcg_csv = pd.read_csv('humanCl.csv')

#create an empty dataframe with the column names in lr_csv
#include gene names that overlap
lr_columnNames = ["Gene Name", "Category", "Description", "Localisation", "Disease Association", "Human ENSEMBL ID"]
lr_df = pd.DataFrame(columns = lr_columnNames)

#gene names in hcg_csv and remove any duplicates
hcg_genes = hcg_csv["Gene Symbol"].drop_duplicates(keep = "first")
#gene names in lr_csv and remove any duplicates
lr_genes = lr_csv["Gene Name"].drop_duplicates(keep = "first")

#gene names that overlap between the two datasets 
sameNames = []
for i in lr_genes:
	for j in hcg_genes:
		if i == j:
			sameNames.append(i)

lr_df["Gene Name"] = sameNames
lr_df["Gene Name"].drop_duplicates(keep = "first")

#set row index based on gene name to get all the info for that row that contains that gene name
lrGeneName_df = lr_csv.set_index("Gene Name")
#cilGeneName_df = cil_df.set_index("Gene Name")

#get the info contained in the "Gene Name", "Description", "Localisation", "Disease Association" and "Human ENSEMBL Id" columns and store it in the appropriate cil_db column
category = []
description = []
localisation = []
diseaseAsso = []
humanEn = []

for gene in lr_df["Gene Name"]:
	category.append(lrGeneName_df.loc[gene]["Category"])
	description.append(lrGeneName_df.loc[gene]["Description"])
	localisation.append(lrGeneName_df.loc[gene]["Localisation"])
	diseaseAsso.append(lrGeneName_df.loc[gene]["Disease Association"])
	humanEn.append(lrGeneName_df.loc[gene]["Human ENSEMBL ID"])

lr_df["Category"] = category 
lr_df["Description"] = description
lr_df["Localisation"] = localisation
lr_df["Disease Association"] = diseaseAsso
lr_df["Human ENSEMBL ID"] = humanEn

#column names for hcg_csv
hcg_colnames = hcg_csv.columns.values

#create dataframe for hcg
hcg_df = pd.DataFrame(columns = hcg_colnames)

#drop the duplicates
hcg_csv["Gene Symbol"] = hcg_csv["Gene Symbol"].drop_duplicates(keep = "first")
hcgGeneName_df = hcg_csv

hcgGeneName_df = hcgGeneName_df.set_index("Gene Symbol")

#for each gene in lr_df add the organism info into hcg_df 
for gene in lr_df["Gene Name"]:
	hcg_df = hcg_df.append(hcgGeneName_df.loc[gene])

lr_df.to_csv("lr_df.csv")
hcg_df.to_csv("hcg_df.csv")
print("Saved lr_df.csv")
print("Saved hcg_df.csv")


