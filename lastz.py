'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 15:52:12
LastEditors: zpliu
LastEditTime: 2022-11-02 20:14:23
@param: 
'''
# *
import os
import pysam
from tempfile import NamedTemporaryFile
from pyliftover import LiftOver 
from gene_coordinate import fasta_seq
lastZ = '/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz'
axtChain = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/axtChain'
chainSwap = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSwap'
chainSort = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSort'


def lastz(seqrchSeq, targetSeq,queryseqFile, targetseqFile, axtFile, chainFile, Chain_changeOrder):
    '''map conserved sequence in subgenome
    '''
    # * write data to tmp file
    queryseqFile.write(seqrchSeq)
    targetseqFile.write(targetSeq)
    queryseqFile.seek(0)
    targetseqFile.seek(0)
    # * align 
    # --filter=identity:85 --filter=coverage:75
    os.system(
        "{} {}[multiple] {}  K=3000 L=3000 H=2000 Y=5000 E=55 T=2 O=600 --filter=identity:80 --filter=coverage:85  --verbosity=10  --format=axt  >{}".format(
            lastZ, targetseqFile.name, queryseqFile.name, axtFile.name
        )
    )

    # * axt to chain
    queryseqFile.seek(0)
    targetseqFile.seek(0)
    axtFile.seek(0)
    os.system(
        "{} -minScore=3000 -linearGap=medium {} -faT -faQ {} {} {}".format(
            axtChain, axtFile.name, targetseqFile.name, queryseqFile.name, chainFile.name
        )
    )
    # * exchange the direction
    chainFile.seek(0)
    os.system(
        "{} {} stdout | {} stdin  {}".format(
            chainSwap, chainFile.name, chainSort, Chain_changeOrder.name
        )
    )
    return None 
def matchRegion(searchId:str,targetId:str,chainTmpFile:NamedTemporaryFile,windowSize=50):
    # ----------------------------------------------------
    # * match the region
    #* searchId: Chr01#100#200#-
    #* targetId: Chr02#100#200#-
    # ----------------------------------------------------
    chainTmpFile.seek(0)
    lo = LiftOver(chainTmpFile.name)
    # Chain_changeOrder.seek(0)
    # with open("test.chain",'w') as File:
    #     File.write(Chain_changeOrder.read())
    # axtFile.seek(0)
    # with open("axtFile.txt",'w') as File:
    #     File.write(axtFile.read())

    # searchId=searchSeq.split("\n")[0].strip(">")
    # targetID=targetSeq.split("\n")[0].strip(">")
    searchStand=searchId.split("#")[-1]
    targetStand=targetId.split("#")[-1]

    
    result = lo.convert_coordinate(searchId, windowSize)
    if targetStand=="+":
        targetStart=int(targetId.split("#")[1])
    else:
        targetStart=int(targetId.split("#")[2])
    targetChrom=targetId.split("#")[0]

    #?映射对应的坐标
    if result==None or len(result)==0:
        return ''
    else:
        if searchStand==targetStand:
            #? same stand
            filterData= [i for i in result if i[2] == "+" and i[0]==targetId]
            matchRegion=[]
            if targetStand=="+":
                for item in filterData:
                    matchRegion.append(
                        "{}#{}#{}".format(
                            targetChrom,targetStart+item[1],targetStand
                        )
                    )
            else:
               for item in filterData:
                matchRegion.append(
                    "{}#{}#{}".format(
                        targetChrom,targetStart-item[1],targetStand
                    )
                )  
        else:
            filterData= [i for i in result if i[2] == "-" and i[0]==targetId]
            #? different stand
            matchRegion=[]
            if targetStand=="+":
                for item in filterData:
                    matchRegion.append(
                        "{}#{}#{}".format(
                            targetChrom,targetStart+item[1],targetStand
                        )
                    )
            else:
                for item in filterData:
                    matchRegion.append(
                        "{}#{}#{}".format(
                            targetChrom,targetStart-item[1],targetStand
                        )
                    ) 
        return ",".join(matchRegion)

    

if __name__ == "__main__":
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    queryseqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    targetseqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    axtFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    chainFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    Chain_changeOrder = NamedTemporaryFile(mode='w+', encoding='utf-8')
    # * 
    searchSeq=fasta_seq(
        genomeObject,Chrom='Ghir_A01',start=186492,end=186692,stand="+"
    )
    targetSeq=fasta_seq(
        genomeObject,Chrom='Ghir_D01',start=163367,end=166050,stand="+"
    )
    lastz(
        searchSeq,targetSeq,queryseqFile=queryseqFile,targetseqFile=targetseqFile,
        axtFile=axtFile,chainFile=chainFile,Chain_changeOrder=Chain_changeOrder,
    )
    queryseqFile.close()
    targetseqFile.close()
    axtFile.close()
    chainFile.close()
    Chain_changeOrder.close()
