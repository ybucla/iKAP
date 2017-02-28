# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:21:03 2017

@author: ybwang
"""
import re,sys,os
from optparse import OptionParser

def main():
    usage = 'usage: python %prog [options] -i MSresultfile -s seqfile -u proteincolumn -p pepcolumn -r ratiocolumn -o output'
    parser = OptionParser(usage)
    parser.add_option('-i', dest='msfile', help='quantification file using MaxQuant or others [Default %default]')
    parser.add_option('-s', dest='seqfile', help='fasta file used for the database search [Default %default]')
    parser.add_option('-o', dest='outdir',default='output', help='output directory [Default %default]')
    parser.add_option('-u', dest='protein',default=0, type='int', help='protein column, 0-based, like "IPI00021812.2" [Default %default]')
    parser.add_option('-p', dest='pep',default=1, type='int', help='pep column, 0-based, like "_HRS(ph)NS(ph)FSDER_" [Default %default]')
    parser.add_option('-r', dest='ratio',default=4, type='int', help='ratio column, 0-based, like "0.38957" [Default %default]')
    (options, args) = parser.parse_args()
    if options.msfile is None or options.seqfile is None or options.protein is None or options.pep is None or options.ratio is None:
        sys.exit("[ERROR] "+parser.get_usage())
    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)
    
    seqdict = readseq(options.seqfile)
    parse(options.msfile,seqdict,options.protein,options.pep,options.ratio)

def parse(infile,seqdict,proteinindex,pepindex,ratioindex, outdir = 'output', rmfirst=True):
    data = {}
    num = {}
    n = 0
    with open(infile,'r') as f:
        for line in f:
            n += 1
            if rmfirst and n == 1: continue
            ele = line.rstrip('\r\n').split("\t")
            id, pep, ratio = ele[proteinindex], ele[pepindex], ele[ratioindex]
            # remove 'K.', '.R', '-', '_', '*', and replace '(ph)' to '#' for convinence
            pep = re.sub(re.compile("(^_)|(_$)"),'',pep)
            pep = re.sub(re.compile("(^\w+\.|^-\.)|(\.\w+$|\.\-$)"),'',pep)
            pep = re.sub(re.compile("@"),'',pep)
            pep = re.sub(re.compile("\(ph\)"),'#',pep)
            pep = re.sub(re.compile("\(\w+\)"),'',pep)
            pep = re.sub(re.compile('\*'),'',pep)
            # if ratio is '' or 'NA', continue
            if ratio == '' or ratio == 'NA': continue
            ratio = float(ratio)
            l = getindex(pep)
            seq = seqdict[id]
            index = seq.find(re.sub(r'#','',pep)) + 1
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

    # output to elm file
    with open(outdir+'/ratio.elm.txt','w') as fout:
        for k in data:
            ratio = data[k] / float(num[k])
            s = k+'\t'+str(1 / ratio)	# attention: ratio reverse only for this paper
            fout.write(s + '\n')
    pdict = {}
    with open(outdir+'/PhosphoPep.txt','w') as fout:
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
	main()
    # seqdict = readseq('seq.txt')
    # parse('data_paper.txt',seqdict,0,1,4)
    # print getPeptideFlank('MSQVQVQVQNPSAALSGSQILNKNQSLLSQPLMSIPSTTSSLPSENAGRPIQNSALPSAS',58,7,7)
