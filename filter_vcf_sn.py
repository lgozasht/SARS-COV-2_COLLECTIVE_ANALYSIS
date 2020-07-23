from ete3 import Tree
import numpy as np
import sys

delList= []
t = Tree(sys.argv[1])
snpDic = {}
sampleDic = {}
finalList = []
mutDic = {}
lineList = []
leafDic = {}
ESPList = []
first = True
#read VCF
with open(sys.argv[2], 'r') as vcf:
    for line in vcf:
        if '##' in line:
            pass
        elif '#' in line:
            sourceList = line.split('\t')
            for e in sourceList:
                    
                ESPList.append(e)
        else:
            lineList.append(line.split('\t'))

#iterate tree
for node in t.traverse("preorder"):
    if node.is_leaf() == True:
        leaf = str(node)[3:]

        leafDic[leaf] = ''
count = 0
finalSourceList = sourceList[0:9]
for i in range(9,len(sourceList)):
    if ESPList[i] in leafDic:
        finalSourceList.append(sourceList[i])


with open("samples_in_latest_tree.txt",'w') as f:
    for source in finalSourceList:
        f.write(source + '\n')





