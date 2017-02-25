import re
def parse(infile,seqdict,proteinindex,pepindex,ratioindex):
    data = {}
    num = {}
    n = 0
    with open(infile,'r') as f:
        for line in f:
            n += 1
            if line.find('Protein') != -1: continue
            ele = line.rstrip('\r\n').split("\t")
            id = ele[0]
            pep = re.sub(re.compile("(^_)|(_$)"),'',ele[1])
            pep = re.sub(re.compile("(^\w+\.|^-\.)|(\.\w+$|\.\-$)"),'',pep)
            pep = re.sub(re.compile("@"),'',pep)
            pep = re.sub(re.compile("\(ph\)"),'#',pep)
            pep = re.sub(re.compile("\(\w+\)"),'',pep)
            pep = re.sub(re.compile('\*'),'',pep)
            if ele[ratioindex] == '': continue
            ratio = float(ele[ratioindex])
            l = getindex(pep)
            seq = seqdict[id]
            index = seq.find(re.sub(r'#','',pep)) + 1
            # print id,pep,l,index
            if len(l) < 1:
                print id,pep,l
            for i in range(len(l)):
                # sites.append(l[i] + index - 1)
                s = l[i] + index - 1
                code = seq[s-1]
                newpep = getPeptideFlank(seq,s,7,7)
                key = id+'\t'+str(s)+'\t'+code+'\t'+newpep
                if key in data:
                    data[key] += ratio
                    num[key] += 1
                else:
                    data[key] = ratio
                    num[key] = 1
                # print id,newpep,s,ratio
    for k in data:
        ratio = data[k] / float(num[k])
        print k+'\t'+str(1 / ratio)
    pdict = {}
    with open('PhosphoPep.txt','w') as fout:
        for k in data:
            pep = k.split('\t')[3]
            center = (len(pep) - 1) / 2
            peptide = pep[0:center] + 'p' + pep[center:]
            peptide = peptide.replace('*','')
            if peptide not in pdict:
                fout.write(peptide+'\n')
                pdict[peptide] = ''

def getPeptideFlank(seq,position,left,right):
    newseq = ('*' * left) + seq + ("*" * right)
    start = left+position - left
    end = right + left + 1 - 1 + start
    return newseq[start-1:end]

def getindex(pep):
    l = []
    for j in range(len(pep)):
        if pep[j] == '#':
            l.append(j-len(l))
    return l

def readseq(seqfile):
    head = ''
    seqdict = {}
    with open(seqfile,'r') as f:
        for line in f:
            if line.find('>') != -1:
                head = line.rstrip().replace('>','')
                seqdict[head] = ''
            else:
                seqdict[head] += line.rstrip()
    return seqdict

if __name__ == '__main__':
    seqdict = readseq('seq.txt')
    parse('data_paper.txt',seqdict,0,1,4)
    # print getPeptideFlank('MSQVQVQVQNPSAALSGSQILNKNQSLLSQPLMSIPSTTSSLPSENAGRPIQNSALPSAS',58,7,7)