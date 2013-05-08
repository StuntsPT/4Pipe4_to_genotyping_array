#!/usr/bin/python2
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#Usage: python2 4Pipe4_to_sequenom.py infile.bam infile.fasta variation_treshold(float: 0-1)

import pysam
import re
from sys import argv

def FASTAtoDict(fasta_file):
    #This will convert the fasta file into a dict like: "name":"seq" and return it
    fasta=open(fasta_file,'r')
    Dict={}
    for lines in fasta:
        if lines.startswith('>'):
            lines=lines.replace('>','')
            name=lines.strip()
            Dict[name]= ''
        else:
            Dict[name] = Dict[name] + lines.upper()
    fasta.close()
    return Dict

def bam_miner(samfile, usable_snps, var_treshold):
    infile = pysam.Samfile(samfile, "rb" )
    var_treshold = float(var_treshold)

    for contig,snps in usable_snps.items():
        for snp in snps:
            snp = int(re.match("(\d+)", snp).group())
            seq_range = set(range(snp - 100, snp) + range(snp + 1, snp + 101))
            for pileupcolumn in infile.pileup(contig):
                if pileupcolumn.pos in seq_range:
                    print "Change"
                    for pileupread in pileupcolumn.pileups:
                        print pileupread.alignment.seq[pileupread.qpos]
                        #TODO: Change print to a counter

        #for ref in infile.references:
            #SNPs = {}
            #raw_snps = ref.split("#")[1:]
            #for i in raw_snps:
                #proc_snps = re.split("(\d+)", i)
                #SNPs[int(proc_snps[1])-1]=set(list(proc_snps[2])) #-1 due to index 0!
                #how_many = len(SNPs)
            #for pileupcolumn in infile.pileup(ref):
                #if pileupcolumn.pos in SNPs:
                    #how_many -= 1
                    #depth = pileupcolumn.n
                    #if depth < min_depth:
                        #Invalidable_snps += 1
                        #continue
                    #bases = set()
                    #for pileupread in pileupcolumn.pileups:
                        #bases.add(pileupread.alignment.seq[pileupread.qpos])
                    #if len(SNPs[pileupcolumn.pos].intersection(bases)) >= 2:
                        #Valid_snps += 1
                    #else:
                        #Invalid_snps += 1
            #if how_many != 0:
                #Invalidable_snps += how_many



def FASTA_miner(fdict):
    usable_snps = {}
    for k,v in fdict.items():
        seq_len = len(v)
        raw_snps = k.split("#")
        contig = raw_snps.pop(0)
        invalid_positions = set(range(100))
        invalid_positions.union(set(range(100-seq_len,seq_len)))
        for i in raw_snps:
            snp_pos = int(re.search("(\d+)", i).group())
            invalid_positions.union(set(range(100-snp_pos,snp_pos+100)))
        for i in raw_snps:
            if int(re.search("(\d+)", i).group()) not in invalid_positions:
                if contig in usable_snps:
                    usable_snps[contig].append(i)
                else:
                    usable_snps[contig] = [i]
    return usable_snps

fdict = FASTAtoDict(argv[2])
usable_snps = FASTA_miner(fdict)
bam_miner(argv[1], usable_snps, argv[3])