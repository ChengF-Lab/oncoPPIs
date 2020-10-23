#!/usr/bin/env python
import re
import sys
import os
import math
from Bio import SwissProt

UniProt_ID_update = {}
for record in SwissProt.parse(open('/Users/junfeizhao/Google Drive/uniprot_sprot_human.dat')):
    for accession in record.accessions:
        UniProt_ID_update[accession] = record.accessions[0]

UniProt_interfaces={}
with open('H_sapiens_interfacesHQ.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        P1 = temp[0]
        P1_IRES = temp[3]
        P2 = temp[1]
        P2_IRES = temp[4]
        if P1 in UniProt_ID_update:
            P1 = UniProt_ID_update[P1]
        if P2 in UniProt_ID_update:
            P2 = UniProt_ID_update[P2]
        if P1 not in UniProt_interfaces:
            UniProt_interfaces[P1] = set()
        if P2 not in UniProt_interfaces:
            UniProt_interfaces[P2] = set()
        if P1_IRES!='[]':
            P1_IRES = P1_IRES[1:-1].split(',')
            for IRES in P1_IRES:
                IRES = IRES.split('-')
                if len(IRES)==1:
                    UniProt_interfaces[P1].add(IRES[0])
                else:
                    UniProt_interfaces[P1].update([str(i) for i in range(int(IRES[0]),int(IRES[1])+1)])
        if P2_IRES!='[]':
            P2_IRES = P2_IRES[1:-1].split(',')
            for IRES in P2_IRES:
                IRES = IRES.split('-')
                if len(IRES)==1:
                    UniProt_interfaces[P2].add(IRES[0])
                else:
                    UniProt_interfaces[P2].update([str(i) for i in range(int(IRES[0]),int(IRES[1])+1)])

#output = open('PPI_sites.txt','w+')
#for Protein in UniProt_interfaces:
#    output.write(Protein + '\t' + ','.join(UniProt_interfaces[Protein]) + '\n')
#output.close()

maf_files = [f for f in os.listdir('/Users/junfeizhao/Google Drive/TCGA-Junfei/Somatic Mutation/') if f.endswith('NonSilent.maf')]
for maf in maf_files:
    Cancer_type = maf.split('.')[1]
    output = open('HQinterfaces_mutations/' + Cancer_type + '_HQinterfaces_mutations.txt',"w+")
    fp = open('/Users/junfeizhao/Google Drive/TCGA-Junfei/Somatic Mutation/' + maf)
    fp.next()
    Header = fp.next()
    #output.write(Header)
    for line in fp:
        temp = line[0:-1].split("\t")
        if temp[8]=='Missense_Mutation':
            gene = temp[0]
            Tumor_ID = temp[15][0:12]
            AA = temp[36]
            AA_pos = re.findall(r'\d+',AA)
            if len(AA_pos) > 0:
                pos = AA_pos[0]
                Protein_ID = temp[67]
                SIFT_score = temp[71]
                PolyPhen_score = temp[72]
                if Protein_ID in UniProt_interfaces and pos in UniProt_interfaces[Protein_ID]:
                    output.write('\t'.join([Tumor_ID,gene,Protein_ID,AA,pos,SIFT_score,PolyPhen_score]) + '\n')
    fp.close()
    output.close()
