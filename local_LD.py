import os
import argparse


parser = argparse.ArgumentParser(description='Identify linked variants in SARS-CoV-2 genomic data') 
parser.add_argument('-v', nargs='?', required=True,
                    help='Path to VCF file')
parser.add_argument('-table', nargs='?', required=True,
                    help='Path to table produced by SARS-COV-2_COLLETIVE_ANALYSIS.py (available at https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter)') 
parser.add_argument('-o', nargs='?', required=True,
                    help='Path to output directory') 
parser.add_argument('-range', nargs='?', required=False, default=10,
                    help='Max range in which LD will be estimated for a given site (Default = 10bp)') 
parser.add_argument('-dependencies', nargs='?', required=False,
                    help="I suggest you download plink2")
#parser.add_argument('-edit_output', nargs='?', required=False,default=False,
#                    help="Add an extra column to final_table.tsv containing ld output")

args = vars(parser.parse_args())

print("Indexing data")
#print(args['edit_output'])
linkageDic = {}
vcfDic = {}
with open(args['v'],'r') as f:
    for line in f:
        if '#' in line:
            pass
        else:
            cols = line.split('\t')
            vcfDic[int(cols[1])] = cols[2]


tableDic = {}
xDic = {}
with open(args['table'],'r') as f1:
    for line in f1:
        if 'Reference' in line:
            header = line
            next(f1)
        else:


            cols = line.strip('\n').split('\t')
            try:
                xDic[int(cols[2])] =  vcfDic[int(cols[2])]
                linkageDic[vcfDic[int(cols[2])]] = {}
                tableDic[str(cols[2])] = cols
            except KeyError:
                print('Warning: position {0} does not exist in provided vcf'.format(str(cols[2])))

for snp in linkageDic:
    if ',' in snp:
        pos = int(snp.split(',')[0][1:-1])
    else:
        pos = int(snp[1:-1])
    for comp in xDic:
        if pos!=comp and (comp > pos-10 and comp < pos+10):
            linkageDic[snp][xDic[comp]] = ' '

finalDic = {}
editDic = {}

print("Performing pairwise comparisons (this might take a while)")
with open('{0}/linked_sites.tsv'.format(args['o']), 'w') as linkFile:
    linkFile.write('snp\tpos\tlinked site|r^2\n')

for snp in linkageDic:
    for comp in linkageDic[snp]:
        os.system("plink2 --allow-extra-chr --vcf {3} --ld {0} {1} > {2}/{0}vs{1}.out".format(snp.replace('U','T'),comp.replace('U','T'),args['o'],args['v']))
        with open("{2}/{0}vs{1}.out".format(snp.replace('U','T'),comp.replace('U','T'),args['o']),'r') as f:
            for line in f:
                if 'r^2' in line:
                    if float(line.strip().split(' ')[2]) > 0.4:
                        if ',' in snp:
                            pos = snp.split(',')[0][1:-1]
                        else:
                            pos = snp[1:-1]
                        if ',' in comp:
                            posComp = comp.split(',')[0][1:-1]
                        else:
                            posComp = comp[1:-1]
                        if snp in finalDic:
                            finalDic[snp][comp]=float(line.strip().split(' ')[2])
                            editDic[pos][posComp]=float(line.strip().split(' ')[2])
                        else:
                            finalDic[snp] = {comp:float(line.strip().split(' ')[2])}
                            editDic[pos] = {posComp:float(line.strip().split(' ')[2])}
        os.system('rm -r {2}/{0}vs{1}.out'.format(snp.replace('U','T'),comp.replace('U','T'),args['o']))
    with open('{0}/linked_sites.tsv'.format(args['o']), 'a') as linkFile:
        if ',' in snp:
            pos = snp.split(',')[0][1:-1]
        else:
            pos = snp[1:-1]


        if snp in finalDic:

            linkFile.write('{0}\t'.format(snp))
            linkFile.write('{0}\t'.format(pos))
            for linkedSite in finalDic[snp]:
                linkFile.write('{0}|{1}\t'.format(linkedSite,str(finalDic[snp][linkedSite])))
            linkFile.write('\n')

        else:
            linkFile.write('{0}\t{1}\tNA\n'.format(snp,pos))
"""
if args['edit_output'] == None:
    with open('{0}/final_table_ld.tsv'.format(args['o']),'w') as f:
        f.write(header)
        for pos in tableDic:
            f.write('\t'.join(tableDic[pos]))
            if pos in editDic:
                f.write('\t')
                for comp in editDic[pos]:
                    f.write('{0}|{1}\t'.format(comp,str(finalDic[pos][comp])))
                f.write('\n')
            else:
                linkFile.write('\tNA\n')


"""
