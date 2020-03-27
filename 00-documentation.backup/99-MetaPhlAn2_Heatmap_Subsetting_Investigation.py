#!/usr/env/bin python

## Stripped down from metaphlan2's metaphlan_hclust_heatmap.py (Jun 2018)

import sys
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from scipy import stats

## Additional investigation
import csv
from pandas import DataFrame

## Set Args

tax_units = "kpcofgs"

fin = "/home/fellows/projects1/calculus_microbiome/evolution/04-analysis/screening/metaphlan2/output/mp2_merged_abundance_table_all_20180628.txt"

xstart = 1
xstop = 9999
ystart = 1
ystop = 9999
percentile = 90 
top = 40
tax_lev = 'g'


## Load data and put into some object (array? Dataframe? Lists?)
mat = [l.strip().split('\t') for l in open( fin ) if l.strip()]

## Extract Genus level
i = tax_units.index(tax_lev) 
mat = [m for i,m in enumerate(mat) if i == 0 or m[0].split('|')[-1][0] == tax_lev or ( len(m[0].split('|')) == i and m[0].split('|')[-1][0].endswith("unclassified"))]

sample_labels=mat[0][xstart:xstop]

m = [(mm[xstart-1],np.array([float(f) for f in mm[xstart:xstop]])) for mm in mat[ystart:ystop]]







## Filter by percentile !!!!
## sorted = sorts list of values
## 		key "A function that would serve as a key or a basis of sort comparison"
## Lambdas are one line functions. They are also known as anonymous functions in some other languages
## 		Lambdas argument: manipulate(argument)
## sort the results of reverse(?)stats.scoreatpercentile(matrix[1],90)

## On each row (taxon), 

m_new = sorted(m,key=lambda x:-stats.scoreatpercentile(x[1],percentile))

## iterate through x == iterate through each entry in M
## apply scoreatpercentile on the iteration, and make minus
## Sort by output of previous step, so 'largest' minus is smallest value so 
## comes first

## For example

for i in range(len(m)):
	-stats.scoreatpercentile(m[i][1], percentile)

## then sort these percentile scores of each taxa indepdent of sample


## Grab the top forty taxa based on the percentile scores
feat_labels = [mm[0].split("|")[-1] for mm in m_new[:top]]

## Filter matrix by these fourty
m_final = [mm[1] for mm in m_new[:top]]
    
D = np.matrix(  np.array( m_final ) )

#return D, feat_labels, sample_labels



## Convert to dataframe with pandas
df = DataFrame.from_records(m)
df.columns = [sample_labels]


## Write files
