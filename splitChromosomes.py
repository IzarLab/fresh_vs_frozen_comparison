inFI = open("/data/Homo_sapiens_assembly19.fasta")
firstline = True
for line in inFI:
    line = line.strip()
    if line.startswith(">"):
        chrname = line.split(' ')[0]
        chrname = chrname[1:len(chrname)]
        if firstline:
            firstline = False
        else:
            outFI.close()
        outFI = open("/data/test_FREEC/"+chrname+".fasta","w")
    outFI.write(line+'\n')
inFI.close()
outFI.close()

inFI = open("/data/Homo_sapiens_assembly19.fasta")
outFI = open("/data/Homo_sapiens_assembly19_good_chrs.fasta",'w')
for line in inFI:
    line = line.strip()
    if line.startswith(">"):
        chrname = line.split(' ')[0]
        chrname = chrname[1:len(chrname)]
    if chrname in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
        outFI.write(line+'\n')
inFI.close()
outFI.close()
