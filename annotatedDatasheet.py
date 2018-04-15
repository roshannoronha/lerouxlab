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

#set row index based on gene name to get all the info for that row that contains that gene name
lrGeneName_df = lr_df.set_index("Gene name")

#get the info contained in the "Gene Name", "Description", "Localisation", "Disease Association" and "Human ENSEMBL Id" columns and store it in the appropriate cil_db column
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

#transfer organisms and ciliopathy values from hcg_df to cil_db
ciliopathyValues = hcg_df.iloc[:, 12:]
for i in range(0, len(orgNames)):
	cil_df[orgNames[i]] = ciliopathyValues.iloc[:,i]

cil_df.to_csv("annotatedTest.csv")
print("Saved annotatedTest.csv")
