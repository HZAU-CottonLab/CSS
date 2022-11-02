'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-10-18 22:09:15
LastEditors: zpliu
LastEditTime: 2022-11-02 11:59:33
@param: 
'''
import pandas as pd
import re
import sys


class Isoform:
    def __init__(self, transcriptId, geneId, start, end, chrom, stand):
        self.id = transcriptId
        self.geneId = geneId
        self.start = start
        self.chrom = chrom
        self.stand = stand
        self.end = end
        self.exons = []

    def set_exon(self, start, end):
        self.exons.append(start)
        self.exons.append(end)

    def get_intron_coordinate(self):
        IntronList = []
        if len(self.exons) == 0:
            return []
        else:
            for i in range(1, len(self.exons)-1, 2):
                IntronList.append(
                    (self.chrom,
                        self.exons[i]+1,
                        self.exons[i+1]-1
                     )
                )
        return IntronList

    def isoform_bed_form(self):
        # * isoform 的bed文件
        return self.chrom, self.id, self.start, self.end, self.geneId, self.stand

    @property
    def info(self):
        return self.chrom, self.start, self.end, self.id, self.stand

    @property
    def gene(self):
        return self.geneId 



def gtf_read_isoform(gtfFile):
    # * isoform dict
    # * PcaBio gtf file
    gtf_annotate = {}
    gtfObject = pd.read_csv(gtfFile, comment="#",
                            header=None, sep="\t", index_col=None)
    for Chrom, itemType, start, end, stand, description in gtfObject[[0, 2, 3, 4, 6, 8]].values:
        if itemType == "transcript" or itemType == 'mRNA':
            isoformId = re.findall(
                'transcript_id \"([^\"]*)\";', description)[0]
            geneID = re.findall('gene_id \"([^\"]*)\";', description)[0]
            gtf_annotate[isoformId] = gtf_annotate.get(isoformId, Isoform(
                isoformId, geneID, start, end, Chrom, stand
            ))
        elif itemType == 'exon':
            gtf_annotate[isoformId].set_exon(start, end)
    return gtf_annotate


class Gene:
    def __init__(self, geneId) -> None:
        self.geneId = geneId
        self.isoformObject = []

    def set_isoform(self, isoformObject):
        self.isoformObject.append(
            isoformObject
        )
    def get_isoformArray(self):
        return self.isoformObject
    @property
    def Isoforms_Id(self):
        return [i.info[3] for i in self.isoformObject]
    @property
    def splic_site(self):
        #* only extract multiple exon isoform
        IsoformArray=[]
        for isoformObject in self.isoformObject:
            chrom, start, end, isoformid, stand=isoformObject.info
            IntonList= isoformObject.get_intron_coordinate()
            if len(IntonList)==0:
                #* exonCount
                #? no intron 
                IsoformArray.append(
                    (isoformid,chrom,start,end,stand,1,';')
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
                    (isoformid,chrom,start,end,stand,len(IntonList)+1,';'.join(spliceArray))
                )
        return IsoformArray


def match_isoformId(description):
    if re.findall('transcript_id \"([^\"]*)\";', description):
        return re.findall('transcript_id \"([^\"]*)\";', description)[0]
    elif re.findall('ID=([^\;]*);?', description):
        matchTxt=re.findall('ID=([^\;]*);?', description)[0]
        #! exon match
        return re.split('\.exon[0-9]+',matchTxt)[0]
    elif re.findall('Parent=([^\;]*);?', description):
        # * exon match
        return re.findall('Parent=([^\;]*);?', description)[0]
    else:
        raise Exception("Isoform Id error")


def match_geneId(description):
    if re.findall('gene_id \"([^\"]*)\";', description):
        return re.findall('gene_id \"([^\"]*)\";', description)[0]
    elif re.findall('Parent=([^\;]*);?', description):
        return re.findall('Parent=([^\;]*);?', description)[0]
    else:
        raise Exception("GeneId error")


def gtf_read_gene(gtfFile):
    # * gene class
    gtf_annotate = {}
    geneDict = {}
    gtfObject = pd.read_csv(gtfFile, comment="#",
                            header=None, sep="\t", index_col=None)
    for Chrom, itemType, start, end, stand, description in gtfObject[[0, 2, 3, 4, 6, 8]].values:
        if itemType == "transcript" or itemType == 'mRNA':
            try:
                geneID = match_geneId(description=description)
                isoformId = match_isoformId(description=description)
            except Exception as error:
                #! raise error
                sys.exit("{} at {}".format(error, description))
            gtf_annotate[isoformId] = gtf_annotate.get(isoformId, Isoform(
                isoformId, geneID, start, end, Chrom, stand
            ))
        elif itemType == 'exon':
            isoformId = match_isoformId(description=description)
            gtf_annotate[isoformId].set_exon(start, end)
    # * merge all isoform
    for isoformId, IsoformObject in gtf_annotate.items():
        geneId = IsoformObject.gene
        geneDict[geneId] = geneDict.get(geneId, Gene(geneId))
        geneDict[geneId].set_isoform(IsoformObject)
    return geneDict
