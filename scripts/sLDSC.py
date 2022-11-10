# ===============================================================
# Step1. Create binary annotation (for ca. 10M snps, split per chromosome) from genomic intervals (Python)   OR from pval of the secondary trait 
# ===============================================================

import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

df_annot_nss = pd.read_table(r'Pruefer2014_MM.txt', header=None, names=['CHR', 'FROM', 'TO'])
df_annot_brain = pd.read_csv(r'H:\Dropbox\analysis\2017_02_February_28_cognition_Neanderthal\3.14\hpabraingenes_noheader.txt', sep='\t')

for chri in range(1, 23):
    df = pd.read_csv(r'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Annot\1000G_Phase3_baselineLD_ldscores\baselineLD.{0}.annot.gz'.format(chri), delim_whitespace=True)
    df = df[['CHR', 'BP', 'SNP', 'CM']].copy()
    
    df_baseline = pd.read_csv(r'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Annot\1000G_EUR_Phase3_baseline\baseline.{0}.annot.gz'.format(chri), delim_whitespace=True)

    t_nss = df_annot_nss[df_annot_nss['CHR'] == 'chr{}'.format(chri)]
    t_nss = IntervalTree.from_tuples(list(zip(t_nss['FROM'], t_nss['TO'])))
    df['NSS'] = [int(bool(t_nss[p])) for p in df['BP']]

    t_brain = df_annot_brain[df_annot_brain['CHR'] == 'chr{}'.format(chri)]
    t_brain = IntervalTree.from_tuples(list(zip(t_brain['FROM'], t_brain['TO'])))
    df['BRAIN'] = [int(bool(t_brain[p])) for p in df['BP']]
    
    df['NSSBRAIN'] = df['NSS'] & df['BRAIN']
    df['CODING'] = df_baseline['Coding_UCSC.bed'] | df_baseline['Intron_UCSC.bed'] | df_baseline['UTR_3_UCSC.bed'] | df_baseline['UTR_5_UCSC.bed']

    print(df.shape, df['NSS'].sum(), df['BRAIN'].sum(), df['NSSBRAIN'].sum(), df['CODING'].sum())
    
    for feature in ['NSS', 'BRAIN', 'NSSBRAIN', 'CODING']:
        df[['CHR', 'BP', 'SNP', 'CM', feature]].to_csv(r'H:\Dropbox\analysis\2018_01_15_NSS\2018_02_26\{0}.{1}.annot.gz'.format(feature, chri), index=False, sep='\t', compression='gzip')


# ===============================================================	
# Step2. Convert binary annotations into LD scores (shell)
# ===============================================================

export LDSC_PY=/mnt/h/GitHub/ldsc/ldsc.py
export LDSCDATA=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR/LDSR_Data
export LDSCANNOT=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR/LDSR_Annot
export NSSANNOT=/mnt/h/Dropbox/analysis/2018_01_15_NSS/2018_02_26
export RESULT=/mnt/h/Dropbox/analysis/2018_01_15_NSS/2018_02_26
export FEATURES='NSS BRAIN NSSBRAIN CODING'
for CHR in {1..22}; do for FEATURE in $FEATURES; do echo "python $LDSC_PY  --l2 --bfile ${LDSCANNOT}/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR --ld-wind-cm 1   --print-snps ${LDSCANNOT}/1000G_Phase3_baselineLD_ldscores/list.txt --annot ${NSSANNOT}/$FEATURE.$CHR.annot.gz --out ${NSSANNOT}/$FEATURE.$CHR"; done; done | parallel -j16

# ===============================================================
# Step3. Run stratified LDSR to compute partitioned heritability (shell)
# ===============================================================

export TASKS='GIANT_BMI_2015_EUR_lift GIANT_HEIGHT_2014_lift UKB_COLLEGE_2016 SSGAC_EDU_2016 CHARGE_COG_2015_lift CTG_INTELLIGENCE_2017  UKB_VNR_2016 UKB_RT_2016 CTG_INTELLIGENCE_2018 PGC_SCZ_2014 PGC_BIP_2016_qc'
for TASK in $TASKS; do for FEATURE in $FEATURES; do echo "python $LDSC_PY --h2 ${LDSCDATA}/${TASK}_noMHC.sumstats.gz --out ${RESULT}/${TASK}.${FEATURE}.partitioned  --ref-ld-chr ${LDSCANNOT}/1000G_EUR_Phase3_baseline/baseline.,${NSSANNOT}/${FEATURE}.    --w-ld-chr ${LDSCANNOT}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.   --overlap-annot    --print-coefficients  --frqfile-chr ${LDSCANNOT}/1000G_Phase3_frq/1000G.EUR.QC.     "; done; done | parallel --dry-run

# ===============================================================
# Step4. Combine resulting files together (python)
# ===============================================================

# Aggregate results together
import glob
import os
import re
import pandas as pd
import numpy as np
import scipy.stats

dir = r'H:\Dropbox\analysis\2018_01_15_NSS\2018_02_26\*.partitioned.results'
files = glob.glob(dir)
df_total = None
for fullfile in files:
    file = os.path.split(fullfile)[1]
    df = pd.read_csv(fullfile, delim_whitespace=True)
    df['file'] = file
    df_total = (df if df_total is None else df_total.append(df))

df_total['TRAIT'] = [x.split('.')[0] for x in df_total['file']]
df_total['FEATURE'] = [x.split('.')[1] for x in df_total['file']]
df_total.to_csv(r'H:\Dropbox\analysis\2018_01_15_NSS\2018_02_26_partitioned.results.csv', index=False, sep='\t')

df_total['Coefficient_p-value'] = scipy.stats.norm.sf(abs(df_total['Coefficient_z-score']))*2  # 2-sided

# ===============================================================
# Step5. Filter table using Excel. (manually)
# ===============================================================
