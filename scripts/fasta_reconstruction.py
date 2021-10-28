#!/usr/bin/python
# -*- coding: utf-8 -*-
#author: Clara Emery
#email: cla.emry@gmail.com 
#last update: August 12th, 2021
import time, argparse, os
from Bio import SeqIO
import utils


def fasta_reconstruct( f_n_fasta_original, f_reconstructed_fasta, threshold=1000): 
#By default if the threshold is not define it will be 1000bp (minium length recognized to carry a gene)
    for record in SeqIO.parse(f_n_fasta_original, 'fasta'): #Biopython parse files as fasta and also genbank etc. more details here: https://biopython.org/wiki/SeqIO
        if len(record.seq) > threshold : #Sorting of all the sequences having a length under the threshold in arg
                SeqIO.write(record, f_plasmids_reconstructed_fasta, 'fasta')

def fasta_reconstruct_plas_meta_sequence( f_n_fasta_original, f_meta_reconstructed_fasta, f_plas_reconstructed_fasta, r_3tools, threshold=1000):
    #Use this function when you have a (read) file with id for sequences identified as a plasmid
    set_3_tools=set()
    for ids_line in r_3tools:
        set_3_tools.add(str(ids_line[0]))  

    for record in SeqIO.parse(f_n_fasta_original, 'fasta'): 
        if len(record.seq) > threshold : #Normally it is already the case, it is just a precaution (and it gets faster too: 14.129sec against 20.095sec )
            if record.id in set_3_tools : 
                SeqIO.write(record, f_plas_reconstructed_fasta, 'fasta')
            # if record.id not in set_3_tools:
            else: 
                SeqIO.write(record, f_meta_reconstructed_fasta, 'fasta')


def fasta_reconstruct_plasmids_ORF_match (f_n_fasta_ORF, f_reconstructed_fasta_ORF, r_alignment, f_parameters):
    d_orf_match=dict()
    for ids_line in r_alignment:
        orf_id=str(ids_line[0])
        if orf_id in d_orf_match:
            d_orf_match[orf_id]+=1
        else:
            d_orf_match[orf_id]=1
    # for k,v in d_orf_match.items():print(k,v)
    nb_record=0
    nb_aligned=0
    nb_not_aligned=0
    for record in SeqIO.parse(f_n_fasta_ORF, 'fasta'): 
        nb_record+=1
        if record.id in d_orf_match : 
            nb_aligned+=1
            for key, value in d_orf_match.items():
                if record.id==key:
                    SeqIO.write(record, f_reconstructed_fasta_ORF, 'fasta') 
                    with open (f_parameters, 'a') as f:
                        f.write(f"{record.id}\t{len(record.seq)}\t{GC(record.seq)}\t{value}\n") 
                        #perhaps round the GC % it will be better to read
        else: nb_not_aligned+=1
    print (nb_record, nb_aligned,nb_not_aligned)

def fasta_reconstruct_plasmids_contigs (f_n_fasta_contig, f_reconstructed_fasta_contig, r_alignment):
    d_contig_match=dict()
    for ids_line in r_alignment:
        contig_id=str(ids_line[0])
        contig_id=contig_id.split("_")
        contig_id=f'{contig_id[0]}_{contig_id[1]}'##############################Verify if it is working or not###################
        if contig_id in d_contig_match:
            d_contig_match[contig_id]+=1
        else:
            d_contig_match[contig_id]=1
    nb_record=0
    nb_aligned=0
    nb_not_aligned=0
    for record in SeqIO.parse(f_n_fasta_contig, 'fasta'): 
        nb_record+=1
        if record.id in d_contig_match : 
            nb_aligned+=1
            for key, value in d_contig_match.items():
                if record.id==key:
                    SeqIO.write(record, f_reconstructed_fasta_contig, 'fasta') 
        else: nb_not_aligned+=1
    print (nb_record, nb_aligned, nb_not_aligned)