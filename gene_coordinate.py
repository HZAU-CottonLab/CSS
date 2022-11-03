'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-11-01 20:34:25
LastEditors: zpliu
LastEditTime: 2022-11-03 14:44:49
@param: 
'''
import pysam
import pandas as pd 
from gtf.isoform_annotate import Gene, gtf_read_gene
import sys 

class Gene_coordinate:
    def __init__(self,geneId,start,end,chrom,stand) -> None:
        self.geneId=geneId
        self.start=start
        self.end=end
        self.chrom=chrom
        self.stand=stand
        self.homoeologous=None
        self.isoformObject=[]
    def set_isoformObject(self,isofromArray):
        self.isoformObject=isofromArray
    def set_homoeologous(self,geneID):
        self.homoeologous=geneID
    @property
    def Isoforms_Id(self):
        return [i.info[3] for i in self.isoformObject]
    @property
    def homoeologousId(self):
        return self.homoeologous
    @property
    def isoform_splic_site(self):
        #* only extract multiple exon isoform
        IsoformArray=[]
        for isoformObject in self.isoformObject:
            chrom, isoformstart, isoformend, isoformid, stand=isoformObject.info
            IntonList= isoformObject.get_intron_coordinate()
            if len(IntonList)==0:
                #* exonCount
                #? no intron 
                IsoformArray.append(
                    (isoformid,chrom,isoformstart,isoformend,stand,1,'')
                )
            else:
                spliceArray=[]
                for Chrom,start,end  in IntonList:
                    spliceArray.append(
                        "{}#{}#{}".format(Chrom,start,stand)
                    ) 
                    spliceArray.append(
                        "{}#{}#{}".format(Chrom,end,stand)
                    )
                IsoformArray.append(
                    (isoformid,chrom,isoformstart,isoformend,stand,len(IntonList)+1,';'.join(spliceArray))
                )
        return IsoformArray
    def gene_sequence(self,genomeObject):
        return  fasta_seq(genomeObject,self.chrom,self.start,self.end,self.stand)

def fasta_seq(fastaObject:pysam.FastaFile,Chrom:str,start:int,end:int,stand:str):
    #* Chr01#1#10#+
    #* Chr01#1#10#-
    seqId='{}#{}#{}#{}'.format(
        Chrom,start,end,stand
    )
    return ">{}\n{}".format(seqId,
        fastaObject.fetch(
        reference=Chrom,start=start-1,end=end
    ) 
    )

#TODO: 获取基因信息以及其转录本信息

def gene_annotate(gtfFile):
    '''
    gtfFile: gtf or gff3
    '''
    geneDict=gtf_read_gene(gtfFile)
    return geneDict

def get_geneObject(gtfFile,geneBedFile,HomoeologFile):
    #-------------------------------------
    #geneBedFile:
    #Ghir_A01        80323913        80324566        +       Ghir_A01G013980
    # Ghir_A01        80322152        80323048        +       Ghir_A01G013970
    # Ghir_A01        60753987        60754357        +       Ghir_A01G013100
    # Ghir_A01        59980270        59983471        +       Ghir_A01G013050
    #HomoeologFile
    # Garb_01G000060  Grai_02G022260
    # Garb_01G000070  Grai_02G022250
    # Garb_01G000080  Grai_02G022240
    # Garb_01G000110  Grai_04G020680
    # Garb_01G000130  Grai_02G022220
    # Garb_01G000140  Grai_02G022210
    # Garb_01G000150  Grai_02G022200
    # Garb_01G000160  Grai_02G022190
    # Garb_01G000170  Grai_02G022180

    #-------------------------------------
    HomoeologDict={}
    Homoeolog=pd.read_csv(HomoeologFile,header=None,index_col=None,sep="\t")
    for geneA,geneB in Homoeolog.values:
        HomoeologDict[geneA]=HomoeologDict.get(geneA,geneB)
        HomoeologDict[geneB]=HomoeologDict.get(geneB,geneA)
    geneIsoformDict=gene_annotate(gtfFile=gtfFile)
    geneBed=pd.read_csv(geneBedFile,header=None,index_col=None,sep="\t")
    GeneDict={}
    #* save the gene info 
    for Chrom,start,end,stand,geneId in geneBed.values:
        geneObject=Gene_coordinate(geneId,start,end,Chrom,stand=stand)
        geneObject.set_homoeologous(
            HomoeologDict.get(geneId,None)
        )
        try:
            geneObject.set_isoformObject(geneIsoformDict.get(geneId).get_isoformArray())
        except AttributeError:
            sys.exit("gene {} is not in gtf".format(geneId))
        GeneDict[geneId] =GeneDict.get(
            geneId,geneObject
        )
    return GeneDict









