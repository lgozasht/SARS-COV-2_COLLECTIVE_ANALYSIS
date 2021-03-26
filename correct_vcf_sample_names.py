vcf = 'filtered_unresolved.vcf.pre'
with open('{0}'.format(vcf.replace('.pre','')),'w') as newFile:
    with open(vcf,'r') as f:
        for line in f:
            if line[0] == '#' and line[1] != '#':
                sp = line.split('\t')
                headers = sp[0:9]
                newLine = '\t'.join([i.split('|')[1] for i in sp if '|' in i])
                newFile.write('{1}\t{0}\n'.format(newLine,'\t'.join(headers)))
            else:
                newFile.write(line)
