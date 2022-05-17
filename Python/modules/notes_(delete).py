# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 22:31:34 2022

@author: yw_ji
"""


    gene_id "ENSG00000137968"; 
    gene_version "12"; 
    transcript_id "ENST00000370859"; 
    transcript_version "3"; 
    exon_number "10"; 
    gene_name "SLC44A5"; 
    gene_source "ensembl_havana"; 
    gene_biotype "protein_coding"; 
    transcript_name "SLC44A5-001"; 
    transcript_source "ensembl_havana"; 
    transcript_biotype "protein_coding"; 
    tag "CCDS"; 
    ccds_id "CCDS44164"; 
    havana_transcript "OTTHUMT00000026823"; 
    havana_transcript_version "3"; 
    protein_id "ENSP00000359896"; 
    protein_version "3"; 
    tag "cds_end_NF"; 
    tag "cds_start_NF"; 
    tag "mRNA_end_NF"; 
    tag "mRNA_start_NF"; 
    tag "basic";
    
    
    gtf = pr.read_gtf(gtfpath, as_df=True) ## takes long time
    gtf = pd.read_table(path, header=None, comment='#', dtype=object)
    
    tmp = gtf[99620:99622].Attribute.str.rstrip(';').str.replace('"','')
    tmp2 = str.extract(r'(?P<name>\w+\s)(?P<value>\w+;\s)')
    tmp.value = tmp.value.str.replace(';','')
    tmp = tmp.reset_index(level=0).pivot(index='level_0',columns='name',values='value')
    gtf2.drop('Attribute', axis=1).join(tmp)
    
    
    col8 = gtf[8].str.split(r'\s"', expand=True)