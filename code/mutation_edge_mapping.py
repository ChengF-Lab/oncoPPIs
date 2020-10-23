#!/usr/bin/env python
import re
import sys
import os
import math
from multiprocessing import Pool
from Bio import SwissProt

UniProt_ID_update = {}
for record in SwissProt.parse(open('/data/Data/UniProt/uniprot_sprot_human.dat')):
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
        interfaces_P1 = set()
        interfaces_P2 = set()
        if P1_IRES!='[]' or P2_IRES!='[]':
            if P1_IRES!='[]':
                P1_IRES = P1_IRES[1:-1].split(',')
                for IRES in P1_IRES:
                    IRES = IRES.split('-')
                    if len(IRES)==1:
                        interfaces_P1.add(IRES[0])
                    else:
                        interfaces_P1.update([str(i) for i in range(int(IRES[0]),int(IRES[1])+1)])
            if P2_IRES!='[]':
                P2_IRES = P2_IRES[1:-1].split(',')
                for IRES in P2_IRES:
                    IRES = IRES.split('-')
                    if len(IRES)==1:
                        interfaces_P2.add(IRES[0])
                    else:
                        interfaces_P2.update([str(i) for i in range(int(IRES[0]),int(IRES[1])+1)])
            UniProt_interfaces[(P1,P2)] = (interfaces_P1,interfaces_P2)

def mutation_mapping(maf):
    Cancer_type = maf.split('.')[1]
    output = open('HQinterfaces_mutation_edges/' + Cancer_type + '_HQinterfaces_mutation_edges.txt',"w+")
    fp = open('/data/Data/TCGA_Somatic_Mutations/Mutect2/' + maf)
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
                mutated_PPI = ['-'.join(PPI) for PPI in UniProt_interfaces if (Protein_ID==PPI[0] and pos in UniProt_interfaces[PPI][0]) or (Protein_ID==PPI[1] and pos in UniProt_interfaces[PPI][1])]
                if len(mutated_PPI)>0:
                    output.write('\t'.join([Tumor_ID,gene,Protein_ID,AA]) + '\t' + ','.join(mutated_PPI) + '\n')
    fp.close()
    output.close()

maf_files = [f for f in os.listdir('/data/Data/TCGA_Somatic_Mutations/Mutect2') if f.endswith('NonSilent.maf')]
p = Pool(8)
p.map(mutation_mapping,maf_files)
