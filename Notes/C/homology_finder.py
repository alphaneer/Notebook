#!/Users/zgxu/anaconda3/envs/python2/bin/python
from __future__ import print_function
import pandas as pd
from pandas import DataFrame, Series
from Bio import SeqIO


gene_name = 'AT3G18780'
gff_path = "/Volumes/Respberry/bioinfo/genome/Athalina/TAIR10/Annotation/TAIR10_GFF3_genes.gff"
genome_path = "/Volumes/Respberry/bioinfo/genome/Athalina/TAIR10/Sequence/TAIR10.fa"

'''
read gff file with pandas's read_table
read reference sequence(fasta format) with biopython
and extract gene sequence from reference based on gff annotation
'''

# read gff
names = ['seqid','source','type','start','end','score','strand','phase','attribute']
gff_file = pd.read_table(gff_path, names = names)


ref = list(SeqIO.parse(genome_path, 'fasta'))
