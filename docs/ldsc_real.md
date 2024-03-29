This usecase describe how to run LDSC analysis (https://github.com/bulik/ldsc) on Morningness and Intelligence summary statistics data. Commands below assume that you've ``git clone`` this github repository (https://github.com/comorment/ldsc) and also https://github.com/comorment/containers repository as described in [Getting started](https://github.com/comorment/containers#getting-started)section. ``$COMORMENT`` environmental variable should points to that parent folder containing cloned repositories.

1.  Export the path of the summary statistics, name this path as ``sumstats_l`` and export path for reference data
```
export sumstats_ld=$COMORMENT/containers/reference/sumstats
export SINGULARITY_BIND="$COMORMENT/ldsc/reference:/REF:ro"
```

2. Copy these uncompressed sumstats to your working directory and Uncompress them if required 
```
cp $sumstats_ld/SavageJansen_2018_intelligence_metaanalysis.txt.gz .
cp $sumstats_ld/Morningness_sumstats_Jansenetal.txt.gz .

gunzip Morningness_sumstats_Jansenetal.txt.gz
gunzip SavageJansen_2018_intelligence_metaanalysis.txt.gz
```

3. Arranging sumstats file for LDSC analysis via  munge_sumstats.py

```
singularity exec --home $PWD:/home $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/munge_sumstats.py \
--sumstats SavageJansen_2018_intelligence_metaanalysis.txt \
--N  2000 \
--out int_munge \
--merge-alleles /REF/w_hm3.snplist

singularity exec --home $PWD:/home $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/munge_sumstats.py \
--sumstats Morningness_sumstats_Jansenetal.txt \
--out mor_munge \
--merge-alleles /REF/w_hm3.snplist \
--signed-sumstats OR,0
```

4. remove .gz extension for munged sumstats

```
mv mor_munge.sumstats.gz mor_munge.sumstats
mv int_munge.sumstats.gz int_munge.sumstats
```

5. Ready to run LDSC analysis

```
singularity exec --home $PWD:/home $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/ldsc.py  \
--rg int_munge.sumstats,mor_munge.sumstats \
--ref-ld-chr /REF/eur_w_ld_chr/ \
--w-ld-chr /REF/eur_w_ld_chr/ \
--out int_mor
```

The succesfull munge_sumstats.py and ldsc.py results shoud look like this:

Munge Intelligence:
![munge1.png](https://raw.githubusercontent.com/comorment/ldsc/main/docs/ldsc_demo/munge1.png)

Munge Morningness:
![munge2.png](https://raw.githubusercontent.com/comorment/ldsc/main/docs/ldsc_demo/munge2.png)

LDSC:
![ldsc.png](https://raw.githubusercontent.com/comorment/ldsc/main/docs/ldsc_demo/ldsc.png)


