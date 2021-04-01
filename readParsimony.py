import time
import os


def parsReader(pars, min_parsimony):
    
    parsimonyDic = {}
    start_time = time.time()
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

                        '''
                        Issue-2
                        '''
                        try:
                            
                            if int(sp[i].split('=')[1]) >= int(min_parsimony):
                                parsimonyDic[pos] = str(sp[i].split('=')[1])


                        except ValueError:
                            pass #print('Warning: possible issue in parsimony output for {0}'.format(sp[0]))
                    except IndexError:
                        pass
                elif int(sp[-1]) >= int(min_parsimony):
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
                        if int(sp[i].split('=')[1]) >= int(min_parsimony):
                            parsimonyDic[pos] = str(sp[i].split('=')[1])
                    except IndexError:
                        pass
                elif int(sp[-1]) >= int(min_parsimony):
                    parsimonyDic[sp[2]] = sp[-1]
    return parsimonyDic
'''
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
                        if int(sp[i].split('=')[1]) >= int(min_parsimony):
                            parsimonyDic[pos] = str(sp[i].split('=')[1])
                    except IndexError:
                        pass
                elif int(sp[-1]) >= int(min_parsimony):
                    parsimonyDic[sp[2]] = sp[-1]


'''
