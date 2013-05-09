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

from __future__ import division
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
    selected_contigs = usable_snps

    for contig,snps in usable_snps.items():
        for snp in snps:
            snp_pos = int(re.match("(\d+)", snp).group())
            seq_range = set(range(snp_pos - 100, snp_pos) + range(snp_pos + 1, snp_pos + 101))
            for pileupcolumn in infile.pileup(contig):
                if pileupcolumn.pos in seq_range:
                    bases = {"A":0,"C":0,"G":0,"T":0}
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.alignment.seq[pileupread.qpos] in bases:
                            bases[pileupread.alignment.seq[pileupread.qpos]] += 1
                    raw_numbers = list(bases.values())
                    most_frequent = max(raw_numbers)
                    var_rate = most_frequent / sum(raw_numbers)
                    if var_rate <= 1 - (var_treshold / 100):
                        if snp in selected_contigs[contig]:
                            selected_contigs[contig].remove(snp)

    #Remove contigs with no SNPs:
    for x in list(selected_contigs.keys()):
        if selected_contigs[x] == []:
            del selected_contigs[x]
    print len(selected_contigs)
    #TODO: Generate the Fastas and we're done!
    #TODO: Something is wrong with the treshold!!!



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