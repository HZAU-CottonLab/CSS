'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-11-01 20:30:23
LastEditors: zpliu
LastEditTime: 2022-11-03 15:37:21
@param: 
'''
from operator import concat
import pysam
from gene_coordinate import fasta_seq, get_geneObject
from lastz import lastz, matchRegion
from tempfile import NamedTemporaryFile
import pandas as pd
import argparse
import numpy as np

def splicingSiteMap(gtfFileA, gtfFileB, geneABedFile, geneBBedFile, genomeFileA, genomeFileB, HomoeologousFile, windowSize=100):
    '''
    gtfFileA: gene annotated file from genomeA
    gtfFileA: gene annotated file from genome B
    genomeFileA: genome sequence File from genome A
    genomeFileB: genome sequence File from genome B
    Homoeologous File:
    # Garb_01G000060  Grai_02G022260
    # Garb_01G000070  Grai_02G022250
    # Garb_01G000080  Grai_02G022240
    # Garb_01G000110  Grai_04G020680
    # Garb_01G000130  Grai_02G022220
    '''
    queryseqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    targetseqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    axtFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    chainFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    Chain_changeOrder = NamedTemporaryFile(mode='w+', encoding='utf-8')

    genomeAObject = pysam.FastaFile(genomeFileA)
    genomeBObject = pysam.FastaFile(genomeFileB)
    Homoeologous = pd.read_csv(
        HomoeologousFile, header=None, index_col=None, sep="\t")
    geneAObjects = get_geneObject(
        gtfFile=gtfFileA, geneBedFile=geneABedFile, HomoeologFile=HomoeologousFile)
    geneBObjects = get_geneObject(
        gtfFile=gtfFileB, geneBedFile=geneBBedFile, HomoeologFile=HomoeologousFile)
    # * only compare betten homoeologous pairs
    #! this can be multi task
    Gene_conservedSite = []
    for geneA, geneB in Homoeologous.values:
    # for geneA, geneB in [('Ghir_A01G000350', 'Ghir_D01G000360')]:
        try:
            singlegeneObjecGeneA = geneAObjects[geneA]
            singlegeneObjecGeneB = geneBObjects[geneB]
        except KeyError:
            # * scaffold gene not in geneBedfile but in Homoeologous File
            # * Skip this gene pairs
            continue
        targeSeqA = singlegeneObjecGeneA.gene_sequence(genomeAObject)
        targeSeqB = singlegeneObjecGeneB.gene_sequence(genomeBObject)
        targetSeq = "{}\n{}".format(targeSeqA, targeSeqB)
        querySeq = ''
        spliceSite_genomeA = []
        spliceSite_genomeB = []
        for isoformData in singlegeneObjecGeneA.isoform_splic_site:
            if isoformData[6]=='':
                #* no intron
                continue
            spliceSite_genomeA += isoformData[6].split(";")
        for isoformData in singlegeneObjecGeneB.isoform_splic_site:
            if isoformData[6]=='':
                continue
            spliceSite_genomeB += isoformData[6].split(";")
        # * extrace_splice sequence
        for spliceData in pd.unique(spliceSite_genomeA):
            Chrom, site, stand = spliceData.split("#")
            querySeq = "{}\n{}".format(
                fasta_seq(genomeAObject, Chrom, int(site) -
                          windowSize, int(site)+windowSize, stand),
                querySeq
            )
        for spliceData in pd.unique(spliceSite_genomeB):
            Chrom, site, stand = spliceData.split("#")
            querySeq = "{}\n{}".format(
                fasta_seq(genomeBObject, Chrom, int(site) -
                          windowSize, int(site)+windowSize, stand),
                querySeq
            )
        if querySeq=='':
            #* no intron sequence
            Gene_conservedSite.append(
                (geneA, geneB,'-','-',0,0,0,'')
            )
            queryseqFile.truncate()
            targetseqFile.truncate()
            axtFile.truncate()
            chainFile.truncate()
            Chain_changeOrder.truncate()
            continue
        lastz(querySeq, targetSeq, queryseqFile=queryseqFile, targetseqFile=targetseqFile,
              axtFile=axtFile, chainFile=chainFile, Chain_changeOrder=Chain_changeOrder)
        # * search the conserved splicing site
        ConservedSplicingSite = []
        for spliceData in pd.unique(spliceSite_genomeA):
            # ? search from geneA
            Chrom, site, stand = spliceData.split("#")
            searchID = '{}#{}#{}#{}'.format(Chrom, int(
                site)-windowSize, int(site)+windowSize, stand)
            targetID = '{}#{}#{}#{}'.format(
                singlegeneObjecGeneB.chrom, singlegeneObjecGeneB.start, singlegeneObjecGeneB.end, singlegeneObjecGeneB.stand)
            MatchRegion = matchRegion(
                searchID, targetID, Chain_changeOrder, windowSize=windowSize)
            if MatchRegion != '':
                # * Chr01:100-Ghir_A01:100,Ghir_A01:200;
                for matchItem in MatchRegion.split(","):
                    ConservedSplicingSite.append(
                        (spliceData, matchItem)
                    )
        # * genome B data
        for spliceData in pd.unique(spliceSite_genomeB):
            # ? search from geneB
            Chrom, site, stand = spliceData.split("#")
            searchID = '{}#{}#{}#{}'.format(Chrom, int(
                site)-windowSize, int(site)+windowSize, stand)
            targetID = '{}#{}#{}#{}'.format(
                singlegeneObjecGeneA.chrom, singlegeneObjecGeneA.start, singlegeneObjecGeneA.end, singlegeneObjecGeneA.stand)
            MatchRegion = matchRegion(
                searchID, targetID, Chain_changeOrder, windowSize=windowSize)
            if MatchRegion != '':
                # * Chr01:100-Ghir_A01:100,Ghir_A01:200;
                for matchItem in MatchRegion.split(","):
                    ConservedSplicingSite.append(
                        (spliceData, matchItem)
                    )
        #! conserved site between genomes A and B
        ConservedSplicingSite = pd.DataFrame(ConservedSplicingSite)
        if ConservedSplicingSite.shape[0]==0:
            Gene_conservedSite.append(
                    (geneA, geneB, "-", "-",
                     "-", "-", 0,"None")
            )
            continue
        # * matchConserved Intron splicing site
        for isoformA in singlegeneObjecGeneA.isoform_splic_site:
            for isoformB in singlegeneObjecGeneB.isoform_splic_site:
                spliceArrayA = isoformA[6].split(";")
                spliceArrayB = isoformB[6].split(";")
                # * conservedSite
                filterData1 = ConservedSplicingSite.loc[(ConservedSplicingSite[0].isin(spliceArrayA)) & (
                    ConservedSplicingSite[1].isin(spliceArrayB))]
                filterData2 = ConservedSplicingSite.loc[(ConservedSplicingSite[1].isin(spliceArrayA)) & (
                    ConservedSplicingSite[0].isin(spliceArrayB))]
                # * geneID isoformId, intronCountA,intronCountB
                # consevrvedSiteCount = filterData1.shape[0]+filterData2.shape[0]
                conservedSiteText=np.array([])
                if filterData1.shape[0]!=0:
                    conservedSiteText=np.hstack(
                        [conservedSiteText,filterData1.apply(lambda x:"{}*{}".format(x[0],x[1]),axis=1).values]
                    )
                if filterData2.shape[0]!=0:
                    conservedSiteText=np.hstack(
                        [conservedSiteText,filterData2.apply(lambda x:"{}*{}".format(x[1],x[0]),axis=1).values]
                    )
            
                #* unique Data
                consevrvedSiteCount=len(np.unique(conservedSiteText))
                conservedSiteText=";".join(np.unique(conservedSiteText))
                if conservedSiteText=="":
                    conservedSiteText='None'
                Gene_conservedSite.append(
                    (geneA, geneB, isoformA[0], isoformB[0],
                     isoformA[5]-1, isoformB[5]-1, consevrvedSiteCount,conservedSiteText)
                )
        queryseqFile.truncate()
        targetseqFile.truncate()
        axtFile.truncate()
        chainFile.truncate()
        Chain_changeOrder.truncate()
    queryseqFile.close()
    targetseqFile.close()
    axtFile.close()
    chainFile.close()
    Chain_changeOrder.close()
    return Gene_conservedSite


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('-g1',help='GTF/GFF3 for genome A')
    parser.add_argument('-g2',help='GTF/GFF3 for genome B')
    parser.add_argument('-b1',help='gene bed file for genome A')
    parser.add_argument('-b2',help='gene bed file for genome A')
    parser.add_argument('-fa1',help='genome sequence file for genome A')
    parser.add_argument('-fa2',help='genome sequence file for genome B')
    parser.add_argument('-orth',help='homoeologous gene pair file')
    parser.add_argument('-w',help='window size, default 100bp flank from splicing site',default=100)
    parser.add_argument('-o',help='output file')
    args=parser.parse_args()
    outArray=splicingSiteMap(
        args.g1,args.g2,args.b1,args.b2,args.fa1,args.fa2,args.orth,int(args.w)
    )
    pd.DataFrame(outArray,columns=[
        'geneId1','geneId2','isoformId1','IsofromId2',
        'IntronCount1','IntronCount2','ConservedSiteCount','conservedSites']).to_csv(args.o,header=True,index=False,sep="\t")
