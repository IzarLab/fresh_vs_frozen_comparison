inFI = open('KAPA.interval_list')
outFI = open('KAPA.interval_list.bed','w')
for line in inFI:
    if line.startswith("chr"):
       words = line.strip().split('\t')
       newline = '\t'.join([words[0][3:], words[1], words[2]])
       outFI.write(newline+'\n')
inFI.close()
outFI.close()
