'''
Date: March 10, 2018
Purpose: Compare ciliated organisms with transition zones to ciliated organisms without transition zones to indentify ciliary genes that are involved witht he transition zone.
"""

import pandas as pd

df = pd.read_excel('TZProteins2.xlsx')

'''
#ciliated with TZ organisms in spreadsheet
#Ciona intestinalis, Danio rerio, Schistosoma mansoni, Strongylocentrotus purpuratus

#ciliated without TZ organisms in spreadsheet
#Plasmodium falciparum, Toxoplasma gondii, Physcomitrella patens

#nonciliated without TZ organisms in spreadsheet
#Saccharomyces cerevisiae, Ustilago maydis, Dictyostelium discoideum, Entamoeba histolytica, Arabidopsis thaliana, Cyanidioschyzon merolae, Phaeodactylum tricornutum, Populus trichocarpa

#store the names of the organisms in the spreadsheet
organism_names = df.columns[12:]
#store all the gene names
gene_names = df['Gene Symbol']

#store ciliated with TZ organisms
cTZ = ["C.intestinalis", "D.rerio", "S.mansoni", "S.purpuratus"]
#store ciliated without TZ organisms

#store nonciliated organisms
ncnTZ = ["S.cerevisiae", "U.maydis", "D.discoideum"]


#new dataframe with gene names as the rows
gene_df = df.set_index('Gene Symbol')

count_ctz = 0
count_ncnTZ = 0

for i in cTZ:
	count_ctz += df.iloc[0, df.columns.get_loc(i)]

for j in ncnTZ:
	count_ncnTZ += df.iloc[0, df.columns.get_loc(j)]

value_ctz = count_ctz/float(len(cTZ))
value_ncnTZ = count_ncnTZ/float(len(ncnTZ))

print("cTZ: " + str(count_ctz) + "/" + str(len(cTZ)) + " = " + str(value_ctz))
print("ncnTZ: " + str(count_ncnTZ) + "/" + str(len(ncnTZ)) + " = " + str(value_ncnTZ))

print(value_ctz + value_ncnTZ)


#names = ["C.intestinalis", "D.rerio"]
#for i in names:
	#print(organism_names.get_loc(i))

'''

#check for ciliated/non-ciliated
#get total number of ciliated/non-ciliated

#get row indices of ciliated with and without TZ
#this is the total number in each catergory
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































