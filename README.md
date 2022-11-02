<!--
 * @Descripttion:
 * @version:
 * @Author: zpliu
 * @Date: 2022-11-01 20:24:52
 * @LastEditors: zpliu
 * @LastEditTime: 2022-11-02 22:02:47
 * @@param:
-->

# Compare splicing site between genomes (CSS)

### required environment

- Python ≥ 3.6
- pandas 
- numpy
- pysam 
- pyliftover

> importance!, **you must to change th**e `lastz.py` **in this package to specify absolute patgs to software like lastz and UCSC et al**. 

```python
#* lastz.py
lastZ = 'software/lastz-1.04.03/lastz-distrib/bin/lastz'
axtChain = 'software/ucsc_kentUtils/v389/bin/axtChain'
chainSwap = 'software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSwap'
chainSort = 'software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSort'

```

### paraments:

+ `-g1` GTF/GFF3 for genome A
+ `-g2` GTF/GFF3 for genome B
+ `b1` gene bed file for genome A
+ `b2` gene bed file for genome B
+ `-fa1` genome sequence file for genome A, **the sequence must be index with** `samtools` 
+ `-fa2` genome sequence file for genome B, **the sequence must be index with** `samtools` 
+ `-orth` homoeologous pairs file 
+ `-w` flank sequence between splicing site, default 100bp
+ `-o` output file

```bash
python main.py  -h 
usage: main.py [-h] [-g1 G1] [-g2 G2] [-b1 B1] [-b2 B2] [-fa1 FA1] [-fa2 FA2] [-orth ORTH] [-w W] [-o O]

optional arguments:
  -h, --help  show this help message and exit
  -g1 G1      GTF/GFF3 for genome A
  -g2 G2      GTF/GFF3 for genome B
  -b1 B1      gene bed file for genome A
  -b2 B2      gene bed file for genome A
  -fa1 FA1    genome sequence file for genome A
  -fa2 FA2    genome sequence file for genome B
  -orth ORTH  homoeologous gene pair file
  -w W        window size, default 100bp flank from splicing site
  -o O        output file

```

header of output 

1. `geneId1` geneID from genome A 
2. `geneId2` geneID from genome B
3. `isoformId1` isoform ID from genome A
4. `isoformId2` isoform ID from genome B
5. `IntronCount1` intron count for isoformId1
6. `IntronCount2` intron count for isoformId2
7. `ConservedSiteCount` count of conserved splicing site between genomes
8. `ConservedSites` conserved splicing sites between genomes





