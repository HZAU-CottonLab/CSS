'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-11-02 21:10:30
LastEditors: zpliu
LastEditTime: 2022-11-02 21:33:44
@param: 
'''
from main import splicingSiteMap 
import pandas as pd 
gtf='/public/home/zpliu/work/Alternative_review/genomeData/Ghirsutum_genome_HAU_v1.1/Ghirsutum_gene_model.gff3'
geneBed='/public/home/zpliu/work/Alternative_review/HomoeologGene/HomoeologGene_v1.1/At_gene.bed'
geneBed1='/public/home/zpliu/work/Alternative_review/HomoeologGene/HomoeologGene_v1.1/Dt_gene.bed'
homoeolog='test/test1.txt'
genomeseq='/public/home/zpliu/work/Alternative_review/genomeData/Ghirsutum_genome_HAU_v1.1/Ghirsutum_genome.fasta' 


geneConservedSite=splicingSiteMap(
    gtf,gtf,geneBed,geneBed1,genomeseq,genomeseq,homoeolog
)
print(pd.DataFrame(geneConservedSite,columns=[
        'geneId1','geneId2','isoformId1','IsofromId2',
        'IntronCount1','IntronCount2','ConservedSiteCount','conservedSites']))