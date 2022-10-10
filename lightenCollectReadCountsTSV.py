inFI = open("/data/test_hg38_rg.sorted.tsv")
linearr = []
for line in inFI:
    line = line.strip()
    linearr.append(line)
inFI.close()

inFI = open("/data/test1_hg38_rg.sorted.tsv")
linearr1 = []
for line in inFI:
    line = line.strip()
    linearr1.append(line)
inFI.close()

if len(linearr)!=len(linearr1):
    nonsense = nonsense+1

outFI = open("/data/test_hg38_rg.sorted.lightened.tsv","w")
for i in range(len(linearr)):
    line = linearr[i]
    line1 = linearr1[i]
    if line[0]!='@' and not line.startswith("CONTIG"):
        words = line.split('\t')
        words1 = line1.split('\t')
        if int(words[3])>500 and int(words1[3])>500:
            outFI.write(line+'\n')
    else:
        outFI.write(line+'\n')
outFI.close()

outFI = open("/data/test1_hg38_rg.sorted.lightened.tsv","w")
for i in range(len(linearr1)):
    line1 = linearr1[i]
    line = linearr[i]
    if line1[0]!='@' and not line1.startswith("CONTIG"):
        words1 = line1.split('\t')
        words = line.split('\t')
        if int(words[3])>500 and int(words1[3])>500:
            #pass
            #print(words1[3])
            outFI.write(line1+'\n')
        # #print(words1[3])
        # if int(words[3])>500 and len(words1[3])>500:
        #     pass
        #     #print("print1")
        #     #outFI.write(line1+'\n')
    else:
        outFI.write(line1+'\n')
outFI.close()

inFI = open('/data/Homo_sapiens_assembly38.preprocessed.interval_list.header')
outFI = open("/data/test_hg38_rg.sorted.lightened.interval_list","w")
for line in inFI:
    line = line.strip()
    outFI.write(line+'\n')
inFI.close()

for i in range(len(linearr1)):
    line1 = linearr1[i]
    line = linearr[i]
    if line1[0]!='@' and not line1.startswith("CONTIG"):
        words1 = line1.split('\t')
        words = line.split('\t')
        if int(words[3])>500 and int(words1[3])>500:
            outFI.write(words[0]+'\t'+words[1]+'\t'+words[2]+'\t'+'+'+'\t'+'-'+'\n')
outFI.close()
