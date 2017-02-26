# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:21:03 2017

@author: ybwang
"""
from optparse import OptionParser
from collections import defaultdict
import numpy as np
from scipy import stats
import sys,os

def main():
    usage = 'usage: python %prog -i ratioElmfile -g iGPSResultfile -o output'
    parser = OptionParser(usage)
    parser.add_option('-i', dest='elmfile', help='quantification ratio elm file [Default %default]')
    parser.add_option('-g', dest='igpsfile', help='result file from iGPS [Default %default]')
    parser.add_option('-o', dest='outdir',default='output', help='output directory [Default %default]')
    (options, args) = parser.parse_args()
    if options.elmfile is None or options.igpsfile is None:
        sys.exit("[ERROR] "+parser.get_usage())
    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)
		
    datalist = parseiGPS(options.elmfile,options.igpsfile)
    ka(datalist,options.outdir)

def parseiGPS(elmfile,igpsfile):
    data = defaultdict(dict)
    # read elm file
    with open(elmfile,'r') as f:
        for line in f:
            ele = line.rstrip().split("\t")
            id,site,code,peptide,ratio = ele
            center = (len(peptide) - 1) / 2
            key1 = peptide[0:center] + 'p' + peptide[center:]
            key2 = '\t'.join(ele)
#            print peptide,center,key1.replace('*','')
            data[key1.replace('*','')][key2] = ''

    # read iGPS result file
    tag = False
    kinaseInfo = {}
    interaction = defaultdict(dict)
    with open(igpsfile,'r') as f:
        for line in f:
            line = line.rstrip()
            if line.find('# ID')  == 0:
                tag = True
            if line.find('# **')  == 0:
                tag = False
            if tag == True and line != '' and line.find('# ID') != 0:
                ele = line.split("\t")
                key1,kid,kname,kfamily = ele[0],ele[6],ele[7],ele[9]
                kinaseInfo[kid] = kname + "\t" + kfamily
                interaction[key1][kid] = ''
    # output
    result = []
    for p in interaction:
        for kid in interaction[p]:
            kinaseinfo = kinaseInfo[kid]
            for l in data[p]:
                result.append(l+'\t'+kid+'\t'+kinaseinfo)
    return result

def ka(datalist,outdir='./'):
    kid_kname = {}
    kinase_treat = defaultdict(float)
    kinase_control = defaultdict(float)
    name = set()
    treat_up_all,control_up_all = 0.0,0.0
    for line in datalist:
        ele = line.split()
        kid,kname = ele[5],ele[6]
        kid_kname[kid] = kname
        ratio = float(ele[4])
        if ratio >= 1:
            ratio = round(ratio)
            kinase_treat[kid] += ratio
            treat_up_all += ratio
        else:
            ratio = round(1 / ratio)
#            ratio = 1 / ratio
            kinase_control[kid] += ratio
            control_up_all += ratio
        name.add(kid)
    # chisquare test
    result = []
    for kid in name:
        ap = kinase_treat[kid] if kid in kinase_treat else 0.0
        an = treat_up_all - ap
        bp = kinase_control[kid] if kid in kinase_control else 0.0
        bn = control_up_all - bp
#        r = scipy.stats.chi2_contingency(np.array([[ap,an],[bp,bn]]))
        ap,an,bp,bn = int(ap),int(an),int(bp),int(bn)
        r = ChiTestSqure(ap,an,bp,bn)
        outline = kid + "\t" + kid_kname.get(kid) + "\t" + str(ap) + "\t" + str(ap / float(ap + an))+ "\t" + str(bp) + "\t" + str(bp / float(bp + bn)) + "\t" + str(r[0]) + "\t" + str(r[1]) + '\t' + str(r[2])
        result.append(outline)
    outfile = outdir+'/ka.txt'
    with open(outfile,'w') as fout:
        fout.write("ProteinId\tName\tNumInGroup1\tPercent\tNumInGroup2\tPercent\tE-ratio\tChiSquare\tPValue\n")
        for r in result:
            fout.write(r+'\n')

def ChiTestSqure(ap,an,bp,bn):
    jz = 0.5
    n1 = ap + an
    n2 = bp + bn
    n3 = ap + bp
    n4 = an + bn
    total = float(ap + an + bp + bn)
    t_ap = n1 * n3 / total
    t_an = n1 * n4 / total
    t_bp = n2 * n3 / total
    t_bn = n2 * n4 / total
    k1 = (abs(ap - t_ap) - jz) * (abs(ap - t_ap) - jz) / float(t_ap)
    k2 = (abs(an - t_an) - jz) * (abs(an - t_an) - jz) / float(t_an)
    k3 = (abs(bp - t_bp) - jz) * (abs(bp - t_bp) - jz) / float(t_bp)
    k4 = (abs(bn - t_bn) - jz) * (abs(bn - t_bn) - jz) / float(t_bn)
    kfValue = k1 + k2 + k3 + k4
    perA = ap / float(n1)
    perB = bp / float(n2)
    er = perA / perB if perB != 0.0 else 'Infinity'
    pvalue = 1-stats.chi2.cdf(kfValue,1)
    result = (er,kfValue,pvalue)
    return result

if __name__ == '__main__':
	main()
    # datalist = parseiGPS('MH.txt','1.iGPS')
    # ka(datalist)
