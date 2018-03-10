"""
Date: March 10, 2018
Purpose: Compare ciliated organisms with transition zones to ciliated organisms without transition zones to indentify ciliary genes that are involved witht he transition zone.
"""

import pandas as pd

df = pd.read_excel('TZProteins2.xlsx')

#get row indices of ciliated with and without TZ
#the values stores represent the total number in each category
cTz_index = df.loc[(df["Ciliated vs Non ciliated"] == "ciliated") & (df["Transition Zone Present"] == "Yes")].index

cnTz_index = df.loc[(df["Ciliated vs Non ciliated"] == "ciliated") & (df["Transition Zone Present"] == "No")].index

#list of gene names
gene_names = list(df)[3:]

#matrix to store probability of a gene being associated with the TZ
tz_matrix = []

#for each gene get value at a index for ciliated/non-ciliated
for i in gene_names:
	cTz_prob = sum(df[i][cTz_index]) / float(len(cTz_index))
	cnTz_prob = sum(df[i][cnTz_index]) / float(len(cnTz_index))
	tz_prob = cTz_prob - cnTz_prob

	tz_matrix.append([i, tz_prob])

#put matrix into a dataframe to look nice
#add column names
tz_matrix = pd.DataFrame(tz_matrix, columns = ["Gene Name", "TZ Probability"])

#sort dataframe by highest to lowest probability
tz_highestprob = tz_matrix.sort_values(by= "TZ Probability", ascending = False)

print(tz_highestprob)

tz_highestprob.to_csv("test.csv")
