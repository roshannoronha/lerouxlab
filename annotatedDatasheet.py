#Written by Roshan Noronha
#April 4, 2018

import pandas as pd
import numpy as np

#read in Leroux/Ritter spreadsheet
lr_df = pd.read_csv('lerouxRitter.csv')

#read in human ciliary genes spreadsheet
hcg_df = pd.read_csv('humanCl.csv')

#get organism names in hcg_df
orgNames = hcg_df.columns.values[12 :,]

#create an empty dataframe with column names
#include gene names that overlap
columnNames = ["Gene Name", "Description", "Localisation", "Disease Association", "Human ENSEMBL Id"]
columnNames.extend(orgNames)
columnNames.extend(["Total Ciliopathy Organisms", "Total Non Ciliopathy Organisms", "Threshold", "Ciliopathy Probability Score", "Ciliopathy Associated"])

cil_df = pd.DataFrame(columns = columnNames)

#get gene names that overlap between the two datasets
sameNames = []
for i in lr_df["Gene name"].values:
	for j in hcg_df["Gene Symbol"]:
		if i == j:
			sameNames.append(i)

#add gene names cil_df
cil_df["Gene Name"] = sameNames

#for each gene name in cil_db get the "Description", "Localisation", "Disease Association" and "Human ENSEMBL Id" from lr_df and add it to the appropriate cil_db column

#set row index based on gene name to get all the info for that row based on the gene name
lrGeneName_df = lr_df.set_index("Gene name")

description = []
localisation = []
diseaseAsso = []
humanEn = []

for i in cil_df["Gene Name"]:
	description.append(lrGeneName_df.loc[i]["Description"])
	localisation.append(lrGeneName_df.loc[i]["Localisation"])
	diseaseAsso.append(lrGeneName_df.loc[i]["Disease Association"])
	humanEn.append(lrGeneName_df.loc[i]["Human ENSEMBL"])

cil_df["Description"] = description
cil_df["Localisation"] = localisation
cil_df["Disease Association"] = diseaseAsso
cil_df["Human ENSEMBL Id"] = humanEn
print(cil_df)	
