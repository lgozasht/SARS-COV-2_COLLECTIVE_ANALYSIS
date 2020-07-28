"""

Identify new possible problematic sites in SARS-CoV-2 genomic data

"""
import os

os.system('rm -r problematic_sites_sarsCov2.vcf')
os.system('wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf')

probDic = {}
newList = []


with open('problematic_sites_sarsCov2.vcf', 'r') as f:
    for line in f:
        if '#' in line:
            pass
        else:
            cols = line.split('\t')
            probDic[cols[1]] = ''


with open('final_table.tsv','r') as f:
    next(f)
    for line in f:
        cols = line.split('\t')

        if cols[12] != 'NA' and int(cols[14])!=0 and int(cols[5])!=0 and int(cols[5]) > 1:
            if float(cols[13])/float(cols[14]) > .5 or float(cols[4])/float(cols[5]) > .5:
                if (cols[9] != 'NA' or cols[16] != 'NA') and str(cols[2]) not in probDic:
                    newList.append(line.strip('\n'))
        elif int(cols[5])!=0 and int(cols[5]) > 1:
            if float(cols[4])/float(cols[5]) > .5:
                if cols[9] != 'NA' and str(cols[2]) not in probDic:
                    newList.append(line.strip('\n'))

with open('new_sites.tsv','w') as f:
    for i in newList:
        f.write(i+'\n')
