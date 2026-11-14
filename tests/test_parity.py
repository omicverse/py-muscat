from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import spearmanr
from pymuscat import aggregate_data,pb_ds
DATA=Path(__file__).parents[1]/'data'
def fixture():
    counts=pd.read_csv(DATA/'counts.csv',index_col=0)
    meta=pd.read_csv(DATA/'metadata.csv',index_col=0)
    a=ad.AnnData(sparse.csr_matrix(counts.T.to_numpy()),obs=meta,var=pd.DataFrame(index=counts.index))
    a.layers['counts']=a.X.copy()
    return a
def test_aggregation():
    a=fixture();pb=aggregate_data(a)
    for key,counts in pb.counts.items():
        expected=pd.read_csv(DATA/f'aggregate_{key}.csv',index_col=0)
        pd.testing.assert_frame_equal(counts,expected,check_dtype=False)
    np.testing.assert_array_equal(pb.n_cells,pd.read_csv(DATA/'n_cells.csv',index_col=0))
    assert pb.n_cells.loc['B','s1']==0
    np.testing.assert_array_equal(a.X.toarray(),a.layers['counts'].toarray())
@pytest.mark.parametrize('method',['edgeR','limma-voom'])
def test_differential_state(method):
    pb=aggregate_data(fixture())
    design=pd.read_csv(DATA/'design.csv',index_col=0)
    contrast=pd.read_csv(DATA/'contrast.csv',index_col=0)
    actual=pb_ds(pb,design,contrast,method=method)['table'].set_index(['cluster_id','gene']).sort_index()
    expected=pd.read_csv(DATA/f'{method}.csv').set_index(['cluster_id','gene']).sort_index()
    assert actual.index.equals(expected.index)
    np.testing.assert_allclose(actual.logFC,expected.logFC,atol=1e-8,rtol=1e-8)
    for col in ['p_val','p_adj.loc','p_adj.glb']:
        assert spearmanr(actual[col],expected[col]).statistic>=.99
        a=set(actual[col].nsmallest(20).index);b=set(expected[col].nsmallest(20).index)
        assert len(a&b)/len(a|b)>=.9
    # Global FDR is across clusters per comparison, not within each cluster.
    assert (actual['p_adj.glb']>=0).all() and (actual['p_adj.glb']<=1).all()
def test_conflicting_sample_labels():
    a=fixture();a.obs.iloc[0,a.obs.columns.get_loc('group_id')]='wrong'
    with pytest.raises(ValueError,match='multiple groups'):aggregate_data(a)

@pytest.mark.parametrize('method',['edgeR','limma-voom'])
def test_public_kang_subset(method):
    counts=pd.read_csv(DATA/'public_counts.csv',index_col=0)
    meta=pd.read_csv(DATA/'public_metadata.csv',index_col=0)
    a=ad.AnnData(sparse.csr_matrix(counts.T.to_numpy()),obs=meta,var=pd.DataFrame(index=counts.index))
    a.layers['counts']=a.X.copy()
    pb=aggregate_data(a)
    d=pd.read_csv(DATA/'public_design.csv',index_col=0)
    ct=pd.read_csv(DATA/'public_contrast.csv',index_col=0)
    actual=pb_ds(pb,d,ct,method=method)['table'].set_index(['cluster_id','gene']).sort_index()
    expected=pd.read_csv(DATA/f'public_{method}.csv').set_index(['cluster_id','gene']).sort_index()
    assert actual.index.equals(expected.index)
    np.testing.assert_allclose(actual.logFC,expected.logFC,atol=1e-8,rtol=1e-8)
    for col in ['p_val','p_adj.loc','p_adj.glb']:
        assert spearmanr(actual[col],expected[col]).statistic>=.99
        x=set(actual[col].nsmallest(20).index);y=set(expected[col].nsmallest(20).index)
        assert len(x&y)/len(x|y)>=.9
