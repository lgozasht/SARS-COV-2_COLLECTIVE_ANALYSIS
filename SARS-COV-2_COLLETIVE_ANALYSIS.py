import os 
import argparse 
import threading
import time 
from threading import Thread 
import copy 
import re


parser = argparse.ArgumentParser(description='Flag lab and country associated variants in SARS-CoV2 Genomes') 
parser.add_argument('-m', nargs='?', required=True,help='Path to GISAID metadata file') 
parser.add_argument('-v', nargs='?', required=True,
                    help='Path to VCF file') 
parser.add_argument('-tree', nargs='?', required=True,
                    help='Path to newick tree') 
parser.add_argument('-o', nargs='?', required=True,
                    help='Path to output directory') 
parser.add_argument('-min_parsimony', nargs='?', required=False, default=4,
                    help='Minimum parsimony for lab and country associations (Default = 4)') 
parser.add_argument('-dependencies', nargs='?', required=False,default=None,
                    help='Check for dependencies')


args = vars(parser.parse_args())



""" Returning an object from a thread """

class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None
    def run(self):
        #print(type(self._target))
        if self._target is not None:
            self._return = self._target(*self._args,
                                        **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return

"""
Define functions
"""

def specificAlleles(pars, vcf, metadata):

    refChrom = 'NC_045512v2'  
    parsimonyDic = {}
    sourceList = []
    sourceDic = {}
    totalParsimoniousmutationCount = 0
    parsimonySourceDic = {}
    aaChangeDic = {}
    globalRefCountDic = {}
    globalAltCountDic = {}
    first = True
    print('Reading parsimony file')
    with open(pars, 'r') as f:
       if pars == 'filtered_unresolved.txt':
            nucs = {'A':' ','G':' ','T':' ','C':' '}
            for line in f:
                sp = line.split('\t')
                if '(hCoV-19' in line[0:10]:
                    pass
                elif 'NC_045512v2' not in line:
                    try:    
                        if ',' in sp[0]:
                            preList = [a for a in sp[0].split(',')[1:] if a in nucs]
                            if sp[0].split(',')[0][-1] in nucs:
                                i = len(preList) + 2
                            else:
                                i = len(preList) + 1

                            pos = sp[0].split(',')[0][1:-1]
                        else:
                            i = 2
                            pos = sp[0][1:-1]
                        
                        try:
                            if int(sp[i].split('=')[1]) >= int(args['min_parsimony']):
                                parsimonyDic[pos] = str(sp[i].split('=')[1])
                        except ValueError:
                            print('Warning: possible issue in parsimony output for {0}'.format(sp[0]))
                    except IndexError:
                        pass
                elif int(sp[-1]) >= int(args['min_parsimony']):
                    parsimonyDic[sp[2]] = sp[-1]
       else:
            for line in f:
                sp = line.split('\t')
                if '(hCoV-19' in line[0:10]:
                    pass
                elif 'NC_045512v2' not in line:
                    try:    
                        if ',' in sp[0]:
                            preList = sp[0].split(',')
                            i = len(preList) + 1
                            pos = preList[0][1:-1]
                        else:
                            i = 2
                            pos = sp[0][1:-1]
                        if int(sp[i].split('=')[1]) >= int(args['min_parsimony']):
                            parsimonyDic[pos] = str(sp[i].split('=')[1])
                    except IndexError:
                        pass
                elif int(sp[-1]) >= int(args['min_parsimony']):
                    parsimonyDic[sp[2]] = sp[-1]



    print('Reading VCF')
    mafDic = {}
    alleleX = {}
    finalDic = {}
    with open(vcf, 'r') as vcf:
        for line in vcf:
            if '##' in line:
                pass
            elif '#' in line:
                sourceList = line.split('\t')
            else:
                cols = line.split('\t')
                (pos, mut) = (cols[1], cols[2])
                if pos in parsimonyDic:
                    total = len(cols[9:])
                    alt = sum(int(x) > 0 for x in cols[9:] if x.isdigit() == True)
                    ref = cols[9:].count('0')
                    #total = int(cols[7].split(';')[1].split('=')[1])
                    preAlt = cols[7].split(';')[0].split('=')[1]
                    mafDic[mut] = {}
                    globalAltCountDic[mut] = {}
                    alleleX[mut] = {}
                    finalDic[mut] = {}
                    #new shit
                    totalAlt = 0
                    if ',' in preAlt:
                        #totalAlt = sum([int(i) for i in preAlt.split(',')])
                        alleles = cols[4].split(',')
                        #preAltCounts = preAlt.split(',')
                        for i in range(0,len(alleles)):
                            finalDic[mut][alleles[i]] = ''
                            globalAltCountDic[mut][alleles[i]] = sum(int(x)== (i+1) for x in cols[9:] if x.isdigit() == True)
                            totalAlt += globalAltCountDic[mut][alleles[i]]
                            alleleX[mut][int(i+1)] = alleles[i]
                    else:
                        globalAltCountDic[mut][cols[4].strip()] = sum(int(x) > 0 for x in cols[9:] if x.isdigit() == True)

                        alleleX[mut][int(1)] = cols[4].strip()
                        finalDic[mut][cols[4].strip()] = ''

                        totalAlt = int(globalAltCountDic[mut][cols[4].strip()])
                    #ref = total - totalAlt
                    globalRefCountDic[mut] = ref
                    for alt in globalAltCountDic[mut]:
                
                        MAF = float(globalAltCountDic[mut][alt])/float(ref + totalAlt)
                        mafDic[mut][alt] = MAF


                    parsimonySourceDic[mut] = {}
                    for i in range(9,len(cols)):
                        epiId = sourceList[i].split('|')[0]
                        if 'EPI_ISL_' in epiId:
                            pass
                        else:
                            epiId = sourceList[i].split('|')[1]
                        if ':' in cols[i]:
                            parsimonySourceDic[mut][epiId] = cols[i].split(':')[0]
                        else:
                            if cols[i].strip() == '.':
                                parsimonySourceDic[mut][epiId]=0
                            else:
                                parsimonySourceDic[mut][epiId] = cols[i]


                    info = cols[7]
                    infoParts = dict([ part.split('=') for part in info.split(';')  ])
                    if ('AACHANGE' in infoParts):
                        aaChangeDic[mut] = infoParts['AACHANGE']

    oriParsCountDic = {}
    subParsCountDic = {}
    subAccessionDic = {}
    oriAccessionDic = {}
    start_time2 = time.time()
    print('Reading metadata')
    accessionToOri = {}
    accessionToSub = {}
    countryAccessionDic = {}

    with open(metadata, 'r') as meta:
        for line in meta:
            sp = line.split('\t')
            if 'strain' == sp[0]:
                indexPosOri = sp.index('originating_lab')
                indexPosSub = sp.index('submitting_lab')
                indexPosCountry = sp.index('country')

            else:
                try:
                    if sp[indexPosOri] not in oriAccessionDic:
                        oriAccessionDic[sp[indexPosOri]] = [sp[2]]
                        accessionToOri[sp[2]] = sp[indexPosOri] 
                    else:
                        oriAccessionDic[sp[indexPosOri]].append(sp[2])
                        accessionToOri[sp[2]] = sp[indexPosOri]


                    if sp[indexPosSub] not in subAccessionDic:
                        subAccessionDic[sp[indexPosSub]] = [sp[2]]
                        accessionToSub[sp[2]] = sp[indexPosOri]

                    else:
                        subAccessionDic[sp[indexPosSub]].append(sp[2])
                        accessionToSub[sp[2]] = sp[indexPosOri]

                    if sp[indexPosCountry] not in countryAccessionDic:
                        countryAccessionDic[sp[indexPosCountry]] = [sp[2]]

                    else:
                        countryAccessionDic[sp[indexPosCountry]].append(sp[2])

                except IndexError:
                    pass
    trashDic = {}
    susDic = {}
    compDic = {}
    otherDic = {}
    altDic = {}
    refDic = {}
    start_time3 = time.time()
    assocLab = {}
    #print(finalDic)
    #printAll = True
    for ori in subAccessionDic:
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = {}
            compDic[snp] = {}
            assocLab[snp] = {}
            trashDic[snp] = {}
            #finalDic[snp] = {}

            for accession in subAccessionDic[ori]:
                try:
                        

                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        #if int(parsimonySourceDic[snp][accession]) >1:
                            #print(int(parsimonySourceDic[snp][accession]))
                        #print(alleleX[snp])
                        if int(parsimonySourceDic[snp][accession]) not in altCount:
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] = 1
                            #print(altCount)
                        else:
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] += 1 
                            #print(altCount)

                except KeyError:
                    #print(altCount)
                    pass
            #print(altCount)

            for allele in altCount:
                #finalDic[snp][allele] = ''


                try:
                    try:
                        if int(globalAltCountDic[snp][allele]) > 0:
                            if (refCount > 0 or altCount[allele] > 0) and (float(altCount[allele])/float(globalAltCountDic[snp][allele]) > .8):
                            #oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refCount,refCount],[globalAltCountDic[snp][allele]-altCount[allele],altCount[allele]]]))
                            #print('This might work')
                                if allele not in compDic[snp] or float(altCount[allele])/float(globalAltCountDic[snp][allele]) > compDic[snp][allele]:
                                    compDic[snp][allele] = float(altCount[allele])/float(globalAltCountDic[snp][allele])
                                    assocLab[snp][allele] = ori
                     
                                    trashDic[snp][allele] = round(float(altCount[allele])/float(globalAltCountDic[snp][allele])*100,2)

                    except ValueError:
                        pass
                except KeyError:
                    pass
        #print(altCount)
    #print(finalDic)
    for ori in oriAccessionDic:
    
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = {}
            for accession in oriAccessionDic[ori]:
                try:
                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        if int(parsimonySourceDic[snp][accession]) not in altCount:
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] = 1 
                        else: 
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] += 1 

                except KeyError:
                    pass
            for allele in altCount:
                #finalDic[snp][allele] = ''

                try:
                    try:
                        if int(globalAltCountDic[snp][allele]) > 0:

                            if (refCount > 0 or altCount[allele] > 0) and (float(altCount[allele])/float(globalAltCountDic[snp][allele])) > .8:
                            #oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refCount,refCount],[globalAltCountDic[snp][allele]-altCount[allele],altCount[allele]]]))

                                if allele not in compDic[snp] or float(altCount[allele])/float(globalAltCountDic[snp][allele]) > compDic[snp][allele]:
                                    compDic[snp][allele] = float(altCount[allele])/float(globalAltCountDic[snp][allele])
                                    assocLab[snp][allele] = ori
                     
                                    trashDic[snp][allele] = round(float(altCount[allele])/float(globalAltCountDic[snp][allele])*100,2)
                                #finalDic[snp] = True
                    except ValueError:
                        pass
                except KeyError:
                    pass
    countryDic = {}
    for ori in countryAccessionDic:
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = {}
            for accession in countryAccessionDic[ori]:
                try:
                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        if int(parsimonySourceDic[snp][accession]) not in altCount:
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] = 1 
                        else: 
                            altCount[alleleX[snp][int(parsimonySourceDic[snp][accession])]] += 1 

                except KeyError:
                    pass
            for allele in altCount:
                try:
                    if altCount[allele] > 0 and int(globalAltCountDic[snp][allele]) > 0:

                        if snp not in countryDic and (float(altCount[allele])/float(globalAltCountDic[snp][allele]))*100.0 > 80.0:
                            countryDic[snp] = {}
                            countryDic[snp][allele] = [ori,str((float(altCount[allele])/float(globalAltCountDic[snp][allele]))*100.0)]
                        elif allele not in countryDic[snp]:
                            countryDic[snp][allele] = [ori,str((float(altCount[allele])/float(globalAltCountDic[snp][allele]))*100.0)]


                except KeyError:
                    pass
    if finalDic is None:
        print('finalDic')
    if countryDic is None:
        print('countryDic')
    if assocLab is None:
        print('assocLab')
    if trashDic is None:
        print('trashDic')
    if globalAltCountDic is None:
        print('globalAltCountDic')
    if mafDic is None:
        print('mafDic')


    return finalDic,countryDic,assocLab,trashDic,globalAltCountDic,mafDic 


def associate(inputType, vcf, pars, meta, missingDic):
    refChrom = 'NC_045512v2'
    parsimonyDic = {}
    sourceList = []
    sourceDic = {}
    dataDic = {}

    totalParsimoniousmutationCount = 0
    parsimonySourceDic = {}
    aaChangeDic = {}
    globalRefCountDic = {}
    globalAltCountDic = {}
    first = True
    print('{0}: Reading parsimony file'.format(inputType))
    with open(pars, 'r') as f:
        for line in f:
            sp = line.split('\t')
            if '(hCoV-19' in line[0:10]:
                pass
            elif 'NC_045512v2' not in line:
                try:    
                    if ',' in sp[0]:
                        preList = sp[0].split(',')  
                        i = len(preList) + 1
                        pos = preList[0][1:-1]
                    else:
                        i = 2
                        pos = sp[0][1:-1]
                    if int(sp[i].split('=')[1]) >= int(args['min_parsimony']):
                        parsimonyDic[pos] = str(sp[i].split('=')[1])
                except IndexError:
                    pass
                    #print(line)
            elif int(sp[-1]) >= int(args['min_parsimony']):
                parsimonyDic[sp[2]] = sp[-1]

    start_time = time.time()
    print('{0}: Reading VCF'.format(inputType))
    mafDic = {}
    if len(missingDic) == 0:
        with open(vcf, 'r') as vcf:
            for line in vcf:
                if '##' in line:
                    pass
                elif '#' in line:
                    sourceList = line.split('\t')
                else:
                    cols = line.split('\t')
                    (pos, mut) = (cols[1], cols[2])
                    if pos in parsimonyDic:
                        #total = int(cols[7].split(';')[1].split('=')[1])
                        total = len(cols[9:])
                        alt = sum(int(x) > 0 for x in cols[9:] if x.isdigit() == True)
                        ref = cols[9:].count('0')

                        globalRefCountDic[mut] = ref
                        globalAltCountDic[mut] = alt
                        MAF = float(alt)/(float(alt)+ float(ref))
                        mafDic[mut] = MAF


                        parsimonySourceDic[mut] = {}
                        for i in range(9,len(cols)):
                            epiId = sourceList[i].split('|')[0]
                            if 'EPI_ISL_' in epiId:
                                pass
                            else:
                                epiId = sourceList[i].split('|')[1]
                            if ':' in cols[i]:
                                parsimonySourceDic[mut][epiId] = cols[i].split(':')[0]
                            else:
                                if cols[i].strip() == '.':
                                    parsimonySourceDic[mut][epiId]=0
                                else:
                                    parsimonySourceDic[mut][epiId] = cols[i]


                        info = cols[7]
                        infoParts = dict([ part.split('=') for part in info.split(';')  ])
                        if ('AACHANGE' in infoParts):
                            aaChangeDic[mut] = infoParts['AACHANGE']
    else:
        with open(vcf, 'r') as vcf:
            for line in vcf:
                if '##' in line:
                    pass
                elif '#' in line:
                    sourceList = line.split('\t')
                else:
                    cols = line.split('\t')
                    (pos, mut) = (cols[1], cols[2])
                    if pos in parsimonyDic:
                        #total = int(cols[7].split(';')[1].split('=')[1])
                        total = len(cols[9:])
                        

                        globalRefCountDic[mut] = 0
                        globalAltCountDic[mut] = 0
                    


                        parsimonySourceDic[mut] = {}
                        for i in range(9,len(cols)):
                            epiId = sourceList[i].split('|')[0]
                            if 'EPI_ISL_' in epiId:
                                pass
                            else:
                                epiId = sourceList[i].split('|')[1]
                            if ':' in cols[i]:
                                parsimonySourceDic[mut][epiId] = cols[i].split(':')[0]
                            else:
                                if cols[i].strip() == '.' or i in missingDic[pos]:
                                    #print(i,cols[i].strip())
                                    #if int(cols[i].strip())>0:
                                    #    test0+=1
                                    parsimonySourceDic[mut][epiId]=0
                                else:

                                    parsimonySourceDic[mut][epiId] = cols[i]
                                    if int(cols[i].strip())>0:
                                        globalAltCountDic[mut]+=1
                                    else:
                                        globalRefCountDic[mut]+=1
                        alt = globalAltCountDic[mut]
                        ref = globalRefCountDic[mut]
                        MAF = float(alt)/(float(alt)+ float(ref))
                        mafDic[mut] = MAF

                        #print('missing alt alleles',test0,'\n','alt alleles',test1)

                        info = cols[7]
                        infoParts = dict([ part.split('=') for part in info.split(';')  ])
                        if ('AACHANGE' in infoParts):
                            aaChangeDic[mut] = infoParts['AACHANGE']

    print('VCF took {0} seconds'.format(time.time() - start_time))

    oriParsCountDic = {}
    subParsCountDic = {}
    subAccessionDic = {}
    oriAccessionDic = {}
    start_time2 = time.time()
    print('Reading metadata')
    accessionToOri = {}
    accessionToSub = {}
    countryAccessionDic = {}

    with open(meta, 'r') as meta:
        for line in meta:
            sp = line.split('\t')
            if 'strain' == sp[0]:
                indexPosOri = sp.index('originating_lab')
                indexPosSub = sp.index('submitting_lab')
                indexPosCountry = sp.index('country')

            else:
                try:
                    if sp[indexPosOri] not in oriAccessionDic:
                        oriAccessionDic[sp[indexPosOri]] = [sp[2]]
                        accessionToOri[sp[2]] = sp[indexPosOri] 
                    else:
                        oriAccessionDic[sp[indexPosOri]].append(sp[2])
                        accessionToOri[sp[2]] = sp[indexPosOri]


                    if sp[indexPosSub] not in subAccessionDic:
                        subAccessionDic[sp[indexPosSub]] = [sp[2]]
                        accessionToSub[sp[2]] = sp[indexPosOri]

                    else:
                        subAccessionDic[sp[indexPosSub]].append(sp[2])
                        accessionToSub[sp[2]] = sp[indexPosOri]

                    if sp[indexPosCountry] not in countryAccessionDic:
                        countryAccessionDic[sp[indexPosCountry]] = [sp[2]]

                    else:
                        countryAccessionDic[sp[indexPosCountry]].append(sp[2])

                except IndexError:
                    pass
    print('metadata took {0} seconds'.format((time.time() - start_time2)))


    print('Working...')
    trashDic = {}
    susDic = {}
    compDic = {}
    otherDic = {}
    altDic = {}
    refDic = {}
    start_time3 = time.time()
    assocLab = {}
    finalDic = {}

    for ori in subAccessionDic:
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = 0
            for accession in subAccessionDic[ori]:
                try:
                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        altCount += 1
                except KeyError:
                   pass
            try:
                try:
                    if int(globalAltCountDic[snp]) > 0:

                        if (refCount > 0 or altCount > 0) and (float(altCount)/float(globalAltCountDic[snp]) > .8):
                       # oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refCount,refCount],[globalAltCountDic[snp]-altCount,altCount]]))
                        #print('This might work')

                            if snp not in compDic or float(altCount)/float(globalAltCountDic[snp]) > compDic[snp]:
                                compDic[snp] = float(altCount)/float(globalAltCountDic[snp])
                                assocLab[snp] = ori
                                trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altCount)/float(globalAltCountDic[snp])*100,2), ori)
                                finalDic[snp] = True

                except ValueError:
                    pass
            except KeyError:
                pass
            if snp not in finalDic:
                finalDic[snp] = True
    for ori in oriAccessionDic:
    
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = 0
            for accession in oriAccessionDic[ori]:
                try:
                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        altCount += 1

                except KeyError:
                    pass
            try:
                try:
                    if int(globalAltCountDic[snp]) > 0:

                        if (refCount > 0 or altCount > 0) and (float(altCount)/float(globalAltCountDic[snp]) > .8):
                
                        #oddsratio, pvalue = stats.fisher_exact(np.array([[globalRefCountDic[snp]-refCount,refCount],[globalAltCountDic[snp]-altCount,altCount]]))
   
                            if snp not in compDic or float(altCount)/float(globalAltCountDic[snp]) > compDic[snp]:
                                compDic[snp] = float(altCount)/float(globalAltCountDic[snp])
                                assocLab[snp] = ori
                                trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altCount)/float(globalAltCountDic[snp])*100,2), ori)
                                finalDic[snp] = True


                except ValueError:
                    pass
            except KeyError:
                pass
            #if printAll == True and snp not in finalDic:
            #    finalDic[snp] = True

    print('Finished {1} association in {0} seconds'.format(time.time() - start_time3,inputType))




    countryDic = {}
    for ori in countryAccessionDic:
        for snp in parsimonySourceDic:
            refCount = 0
            altCount = 0
            for accession in countryAccessionDic[ori]:
                try:
                    if int(parsimonySourceDic[snp][accession]) == 0:
                        refCount += 1
                    elif int(parsimonySourceDic[snp][accession]) > 0:
                        altCount += 1
                except KeyError:
                    pass
            try:
                if altCount > 0 and int(globalAltCountDic[snp]) > 0:

                    if snp not in countryDic and (float(altCount)/float(globalAltCountDic[snp]))*100.0 >= 80.0:
                        countryDic[snp] = [ori,str((float(altCount)/float(globalAltCountDic[snp]))*100.0)]
                        if snp not in trashDic:
                            trashDic[snp] = '{0}% of alternate allele calls stem from {1}'.format(round(float(altDic[ori][snp])/float(globalAltCountDic[snp])*100,2), ori)

            except KeyError:
                pass




    primerOverlap = {}

    primerTrackList = []
    primerDic = {}

    with open('primers.txt', 'r') as primeFile:


        for line in primeFile:
 
            for snp in finalDic:
  
                sp = line.split('\t')
                if ',' in snp:
                    pos = snp.split(',')[0][1:-1]
                else:
                    pos = snp[1:-1] 
                if int(pos) > int(sp[1]) and  int(pos) < int(sp[2]):
                    if snp not in primerDic:

                        primerOverlap[snp] = sp[3]
                    else:
                        primerDic[snp] += (',' + sp[3])


                elif int(pos) > (int(sp[1]) - 10) and int(pos) < (int(sp[2]) + 10):
                    primerTrackList.append(line)
                    if snp not in primerDic:
                        primerDic[snp] = sp[3]
                    else:
                        primerDic[snp] += (',' + sp[3])
    for snp in finalDic:
        primer = primerDic[snp] if snp in primerDic else 'NA'
        primeroverlap = primerOverlap[snp] if snp in primerOverlap else 'NA'
   
        countryOri = countryDic[snp][0] if snp in countryDic else 'NA'
        countryAt = countryDic[snp][1] if snp in countryDic else 'NA'

        if ',' in snp:
            pos = snp.split(',')[0][1:-1]
        else:
            pos = snp[1:-1]
        ori = countryDic[snp] if snp in countryDic else 'NA'
        assoc = assocLab[snp] if snp in assocLab else 'NA'
        assocValue = trashDic[snp].split(' ')[0] if snp in assocLab else 'NA' 
        dataDic[snp] = {'refChrom':refChrom, 'start':str(int(pos)-1), 'stop':pos,'snp':snp.replace('T', 'U'),
                                'parsimonyDic':parsimonyDic[pos].strip('\n'), 'globalAltCountDic':str(globalAltCountDic[snp]),
                                'mafDic':str(mafDic[snp]),'primeroverlap':primeroverlap,'primer':primer,'countryOri':countryOri,'countryAt':countryAt,
                                'assoc':assoc,'assocValue':assocValue}
    print(dataDic)
    return dataDic, finalDic



def replaceGenotype(cols,replacements,alts):
    try:
        #replacements['.'] = '0'
        #print(replacements)
        if len(replacements) > 0:
            genotypes = re.sub('({})'.format('|'.join(map(re.escape, replacements.keys()))), lambda m: replacements[m.group()], ','.join(cols[9:]))
        else:
            return cols

    except KeyError:
        print(altDic,cols[0:9])
    return [str(x) for x in cols[0:9]+genotypes.split(',')]

def replaceGenotypeSingle(cols):
    
    genotypes = ['1' if (str(x)!='0' and str(x)!='.') else x for x in cols[9:len(cols)]]
    
    return [str(x) for x in cols[0:9]+genotypes]


def resolve_alt():
    print('Resolving with alternates')

    if os.path.isfile('filtered_resolved_alt.txt') == True and os.stat("filtered_resolved_alt.txt").st_size > 0:
        print('File exists... Moving on')
    else:

        ambDic = {"N": {"A":'', "C":'', "G":'', "T":''},
         "R": {"A":'', "G":''},
         "Y": {"T":'', "C":''},
         "K": {"G":'', "T":''},
         "M": {"A":'', "C":''},
         "W": {"A":'', "T":''},
         "S": {"G":'', "C":''},
         "B": {"C":'', "G":'', "T":''},
         "D": {"A":'', "G":'', "T":''},
         "H": {"A":'', "C":'', "T":''},
         "V": {"A":'', "C":'', "G":''}}
        normalDic = {"A":'', "C":'', "G":'', "T":'' }
        with open('filtered_resolved_alt.vcf', 'w') as resolved:

            with open('filtered_unresolved.vcf', 'r') as f:
                for line in f:
                    if '##' in line:
                        resolved.write(line)

                    elif '#' in line:
                         resolved.write(line)
                    else:
                        alts = {}
                        cols = line.split('\t')
                        altList = cols[4].split(',')
                        snpList = cols[2].split(',')
                        altOrder = []
                        snpOrder = []
                        altReplaceDic = {}
                        reorganizeAlleleDic = {}
                        altCounts = cols[7].split(';')[0].split('=')[-1].split(',')
                        ref = 'AN={0}'.format(str(len(cols[9:])))
                        foundAlt = False
                        found = False
                        for i in range(0,len(altList)):
                            if altList[i] in normalDic:
                                
                                altOrder.append(altList[i])
                                snpOrder.append(snpList[i])
                                alts[altList[i]] = int(altCounts[i])
                                reorganizeAlleleDic[str(i+1)] = str(altOrder.index(altList[i])+1)
                                foundAlt = True
                        updatedCols = cols
                        if len(reorganizeAlleleDic.keys()) == len(altList):
                            resolved.write(line)
                            continue
                        elif foundAlt == True:
                    
                            for i in range(0,len(altList)):
                                if altList[i] not in normalDic:
                                    found = False

                                    for k in range(0,len(altList)):
                                        if altList[k] in ambDic[altList[i]] and altList[k] in normalDic:
   
                                            alts[str(altList[k])] += int(altCounts[i])
                                            altReplaceDic[str(i+1)] = str(k+1)

                                            found = True
                                            break

                                    if found==False:
                                        altReplaceDic[str(i+1)] = str(0)
                                    

                            updatedCols = replaceGenotype(updatedCols,altReplaceDic,alts)
                            if len(reorganizeAlleleDic) > 0:
                                updatedCols = replaceGenotype(updatedCols,reorganizeAlleleDic,alts)
                            else:
                                updatedCols = cols
                            error = False
                            d = 0
                            for i in range(0,len(altOrder)):
                                if error == False:
                                    v = i - d
                                    try:
                                        L = updatedCols.index(str(v+1))
                                    except ValueError:
                                        
                                        del alts[altOrder[v]]
                                        altOrder.remove(altOrder[v])
                                        errorFix = {}
                                        error = True
                                        for k in range(v,len(altOrder)):
                                            errorFix[str(k+2)]=str(k+1)
                                        d+=1

                            if error == True:
                                updatedCols = replaceGenotype(updatedCols,errorFix,alts)

                            
                            
                            newCols = updatedCols

                        else:
                            for alt in ambDic[altList[0]]:
                                if alt not in cols[3]:

                                    altCount = sum([int(x) for x in altCounts])
                                    altFinal = alt
                                    snpFinal = '{0}{1}'.format(cols[2][0:-1],alt)
                                    newCols = replaceGenotypeSingle(cols)

                                    break

                        if foundAlt == True:                     
                            newLine = '\t'.join(newCols).replace(newCols[4],','.join(altOrder))
                            newLineCorrectedCounts = newLine.replace(newCols[7],'AC={0};{1}'.format(','.join([str(alts[x]) for x in altOrder]),str(ref)))
                            newLineCorrectedSnps = newLineCorrectedCounts.replace(newCols[2],','.join(snpOrder))
                        else:
                         
                            newLine = '\t'.join(newCols).replace(newCols[4],altFinal)
                            newLineCorrectedCounts = newLine.replace(newCols[7],'AC={0};{1}'.format(str(altCount),str(ref)))
                            newLineCorrectedSnps = newLineCorrectedCounts.replace(newCols[2],snpFinal)+'\n'
                        try:
                            if len(altFinal) == 0 or altFinal == ' ' or len(altOrder)==0:
                                pass
                            else:
                                resolved.write(newLineCorrectedSnps)
                        except UnboundLocalError:

                            resolved.write(newLineCorrectedSnps)


        os.system('./strain_phylogenetics/build/find_parsimonious_assignments --tree {0} --vcf filtered_resolved_alt.vcf > filtered_resolved_alt.txt'.format(args['tree']))


def resolve_ref():
    print('Resolving with reference')
    if os.path.isfile('filtered_resolved_ref.txt') == True and os.stat("filtered_resolved_ref.txt").st_size > 0:
        print('File exists... Moving on')
    else:
        os.system('./strain_phylogenetics/build/find_parsimonious_assignments --tree {0} --vcf filtered_unresolved.vcf --print-vcf > filtered_resolved_ref.vcf'.format(args['tree']))
        os.system('./strain_phylogenetics/build/find_parsimonious_assignments --tree {0} --vcf filtered_resolved_ref.vcf > filtered_resolved_ref.txt'.format(args['tree']))

def unresolved():
    print('Computing parsimony from unresolved vcf')
    if os.path.isfile('filtered_unresolved.txt') == True and os.stat("filtered_unresolved.txt").st_size > 0:
        print('File exists... Moving on')
    else:

        os.system('./strain_phylogenetics/build/find_parsimonious_assignments --tree {0} --vcf filtered_unresolved.vcf > filtered_unresolved.txt'.format(args['tree']))

def missingSamples(vcf):
    missingDic = {}
    with open(vcf, 'r') as vcf:
        for line in vcf:
            if '#' in line:
                pass
            else:
                cols = line.split('\t')
                pos = cols[1]
                missingDic[pos] = get_index_positions(cols[9:],'.')
    return missingDic

def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos+9)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

if args['dependencies'] != None:
    print('Checking for parsimonious assignment software')
    if os.path.isdir('./strain_phylogenetics') == False:
        print('"strain_phylogenetics" must be in your current working directory. You can download it at https://github.com/yatisht/strain_phylogenetics') 
    if os.path.isfile('remove_samples.pl') == False:
        print('./remove_samples.pl must be in your current working directory. You can download it at https://github.com/lgozasht/COVID-19-Lab-Specific-Bias-Filter')

print('Identifying included samples')
if os.path.isfile('samples_in_latest_tree.txt') == True and os.stat("samples_in_latest_tree.txt").st_size > 0:
    print('File exists... Moving on')
else:
    os.system('python3.6 filter_vcf_sn.py {1} {0}'.format(args['v'], args['tree']))

print('Filtering samples from VCF that do not exist in the provided newick tree')

if os.path.isfile('filtered_unresolved.vcf') == True and os.stat("filtered_unresolved.vcf").st_size > 0:
    print('File exists... Moving on')
else:
    os.system("perl remove_samples.pl {0} > filtered_unresolved.vcf".format(args['v']))

print("Indexing missing data")
missingDic = missingSamples('filtered_unresolved.vcf')

print('Resolving ambiguities')

    
t1 = threading.Thread(target=resolve_ref) 
t2 = threading.Thread(target=resolve_alt)
t3 = threading.Thread(target=unresolved) 


t1.start() 
t2.start() 
t3.start()

t1.join()
t2.join()
t3.join() 
    
print('Performing associations')

assocAlt = ThreadWithReturnValue(target = associate, args = ('alternate_resolved','filtered_resolved_alt.vcf',
                                                                 'filtered_resolved_alt.txt',args['m'], {}))
assocRef = ThreadWithReturnValue(target = associate, args = ('reference_resolved','filtered_resolved_ref.vcf',
                                                                 'filtered_resolved_ref.txt',args['m'],missingDic))
alleleSpcecific = ThreadWithReturnValue(target = specificAlleles, args = ('filtered_unresolved.txt','filtered_unresolved.vcf',args['m']))

assocAlt.start()
assocRef.start()
alleleSpcecific.start()

altDic,finalDicAlt = assocAlt.join()
refDic,finalDicRef = assocRef.join()
finalAlleles,countryDic,assocLab,trashDic,globalAltCountDic,mafDic = alleleSpcecific.join()

refStarts = {}
finalAlleles2 = {} 
countryDic2= {}
assocLab2 = {} 
trashDic2={}
globalAltCountDic2 = {}
mafDic2 = {}
refDic2 = {}

for snp in finalDicRef:
    if ',' in snp:
        pos = str(snp.split(',')[0][1:-1])
    else:
        pos = str(snp[1:-1])
    try:
        refStarts[pos] = finalDicRef[snp]
        refDic2[pos] = refDic[snp]
    except KeyError:
        refStarts[pos] = finalDicRef[snp.split(',')[0]]
        refDic2[pos] = refDic[snp]


for snp in finalAlleles:
    if ',' in snp:
        pos = str(snp.split(',')[0][1:-1])
    else:
        pos = str(snp[1:-1])

    finalAlleles2[pos] = finalAlleles[snp] 
    countryDic2[pos]= countryDic[snp] if snp in countryDic else 'NA'
    assocLab2[pos] = assocLab[snp] if snp in assocLab else 'NA'
    trashDic2[pos] = trashDic[snp]  if snp in trashDic else 'NA'
    globalAltCountDic2[pos]=  globalAltCountDic[snp] if snp in globalAltCountDic else 'NA'
    mafDic2[pos] = mafDic[snp] if snp in mafDic else 'NA'
    

with open('{0}/final_table.tsv'.format(args['o']), 'w') as f:
    f.write('\t'.join(['Reference','Start','Stop','Snp',
                       'Alt Resolved Parsimony','MAC','MAF','Primer Overlap',
                       'Primer Vicinity','Country','Country Association','Lab','Lab Association',
                       'Ref Resolved Parsimony','MAC','MAF',
                       'Country','Country Association','Lab','Lab Association','\n']))
    for snp in finalDicAlt:
        if ',' in snp:
            pos = str(snp.split(',')[0][1:-1])
        else:
            pos = str(snp[1:-1])

        if pos in refStarts:


            f.write('\t'.join([altDic[snp]['refChrom'],altDic[snp]['start'], altDic[snp]['stop'], altDic[snp]['snp'],
                           altDic[snp]['parsimonyDic'], altDic[snp]['globalAltCountDic'],
                           altDic[snp]['mafDic'], altDic[snp]['primeroverlap'],
                           altDic[snp]['primer'], altDic[snp]['countryOri'], altDic[snp]['countryAt'], altDic[snp]['assoc'], altDic[snp]['assocValue'], 
                           refDic2[pos]['parsimonyDic'], refDic2[pos]['globalAltCountDic'],
                           refDic2[pos]['mafDic'],
                           refDic2[pos]['countryOri'], refDic2[pos]['countryAt'], refDic2[pos]['assoc'], refDic2[pos]['assocValue']]))
        else:

            f.write('\t'.join([altDic[snp]['refChrom'],altDic[snp]['start'], altDic[snp]['stop'], altDic[snp]['snp'],
                           altDic[snp]['parsimonyDic'], altDic[snp]['globalAltCountDic'],
                           altDic[snp]['mafDic'], altDic[snp]['primeroverlap'],
                           altDic[snp]['primer'], altDic[snp]['countryOri'], altDic[snp]['countryAt'], altDic[snp]['assoc'], altDic[snp]['assocValue'], 
                           'NA', 'NA',
                           'NA', 'NA',
                           'NA', 'NA', 'NA']))

        if pos in finalAlleles2:
            f.write('\t')
            for allele in finalAlleles2[pos]:
                if pos in countryDic2:
                    try:
                        countryOri = countryDic2[pos][allele][0] if allele in countryDic2[pos] else 'NA'
                        countryAt = countryDic2[pos][allele][1] if allele in countryDic2[pos] else 'NA'
                    except TypeError:
                        pass
                else:
                    countryOri = 'NA'
                    countryAt = 'NA'
                if pos in assocLab2:
                    assoc = assocLab2[pos][allele] if allele in assocLab2[pos] else 'NA'
                else:
                    assoc = 'NA'
                if pos in trashDic2:
                    assocValue = str(trashDic2[pos][allele]) if allele in trashDic2[pos] else 'NA'
                else:
                    assocValue = 'NA' 
                f.write('|'.join([allele,str(globalAltCountDic2[pos][allele]),str(mafDic2[pos][allele]), countryOri, str(countryAt),str(assoc),str(assocValue)]) + '\t')
        f.write('\n')

