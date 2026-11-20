python - <<'PY'
from pathlib import Path
import pandas as pd, numpy as np
oldp=Path('/data/projects/punim1236/liver_scRNAseq/AAA_complete_altogether.csv')
newp=Path('/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params/03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv')
old=pd.read_csv(oldp)
old=old[old.source.astype(str).str.strip().isin(['cytotoxic CD4','naive CD4']) & (old.target=='MNP infected')].copy()
old['source_clean']=old.source.str.strip()
# one row per pair/source; retain highest prob if duplicate
old=old.sort_values('prob').drop_duplicates(['interaction_name','source_clean'],keep='last')
wide=old.pivot(index='interaction_name',columns='source_clean',values='prob').reindex(columns=['cytotoxic CD4','naive CD4']).fillna(0)
wide['old_log_cyto']=np.log1p(wide['cytotoxic CD4'])
wide['old_log_naive']=np.log1p(wide['naive CD4'])
wide['old_signed']=np.where(wide['cytotoxic CD4']>=wide['naive CD4'],wide['old_log_cyto'],-wide['old_log_naive'])
new=pd.read_csv(newp).set_index('interaction_name')
comp=wide.join(new[['prob_infected','prob_naive','log_prob_infected','log_prob_naive','signed_log_prob','higher_condition']],how='outer')
for c in ['prob_infected','prob_naive','log_prob_infected','log_prob_naive','signed_log_prob']:
    comp[c]=comp[c].fillna(0)
comp['new_signed']=comp['signed_log_prob']
comp['old_side']=np.where(comp.old_signed>0,'infected/cytotoxic',np.where(comp.old_signed<0,'naive','absent'))
comp['new_side']=np.where(comp.new_signed>0,'infected/cytotoxic',np.where(comp.new_signed<0,'naive','absent'))
comp['delta_signed']=comp.new_signed-comp.old_signed
oldset=set(wide.index); newset=set(new.index)
print('old_unique_pairs',len(oldset))
print('new_unique_pairs',len(newset))
print('overlap',len(oldset&newset))
print('old_only',len(oldset-newset))
print('new_only',len(newset-oldset))
over=comp.loc[sorted(oldset&newset)].copy()
agree=(over.old_side==over.new_side).sum()
print('direction_agreement_overlap',agree,'/',len(over),f'({agree/len(over):.1%})')
print('spearman_signed_overlap',over.old_signed.corr(over.new_signed,method='spearman'))
print('pearson_signed_overlap',over.old_signed.corr(over.new_signed,method='pearson'))
print('median_abs_delta_signed_overlap',over.delta_signed.abs().median())
print('\nTop old stronger infected/cytotoxic:')
print(over.sort_values('old_signed',ascending=False)[['old_signed','new_signed','delta_signed']].head(10).to_string())
print('\nTop old stronger naive:')
print(over.sort_values('old_signed')[['old_signed','new_signed','delta_signed']].head(10).to_string())
print('\nLargest absolute changes among overlap:')
print(over.assign(absdelta=over.delta_signed.abs()).sort_values('absdelta',ascending=False)[['old_signed','new_signed','delta_signed']].head(15).to_string())
print('\nOld-only notable:')
print(comp.loc[sorted(oldset-newset)].assign(absold=lambda x:x.old_signed.abs()).sort_values('absold',ascending=False)[['old_signed','new_signed']].head(15).to_string())
print('\nNew-only notable:')
print(comp.loc[sorted(newset-oldset)].assign(absnew=lambda x:x.new_signed.abs()).sort_values('absnew',ascending=False)[['old_signed','new_signed']].head(15).to_string())
# save comparison
out=Path('/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params/MNP2_vs_original_combinedMNP_comparison.csv')
comp.reset_index(names='interaction_name').to_csv(out,index=False)
print('saved',out)
PY