"""muscat pseudobulk subset, using existing pyedger/python-limma engines."""
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class Pseudobulk:
    counts: dict
    n_cells: pd.DataFrame
    samples: pd.DataFrame


def aggregate_data(adata, *, cluster_key='cluster_id', sample_key='sample_id',
                   group_key='group_id', layer='counts'):
    """Sum counts per cluster x biological sample, retaining absent combinations.

    Each sample must have exactly one group. Returns one genes-by-samples table
    per cluster, cell counts and sample metadata; never modifies AnnData.
    This implements aggregateData(fun='sum', scale=FALSE, by=two keys).
    """
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('Cell and gene names must be unique.')
    for key in (cluster_key,sample_key,group_key):
        if key not in adata.obs or adata.obs[key].isna().any():raise ValueError(f'Missing {key} labels.')
    if layer is not None and layer not in adata.layers:raise ValueError('Missing count layer.')
    x=adata.X if layer is None else adata.layers[layer]
    values=x.data if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all() or np.any(values<0) or not np.allclose(values,np.rint(values),rtol=0,atol=1e-8):
        raise ValueError('Aggregation requires finite nonnegative integer counts.')
    samples=sorted(pd.unique(adata.obs[sample_key]));clusters=sorted(pd.unique(adata.obs[cluster_key]))
    if not samples or not clusters:raise ValueError('No samples/clusters remain.')
    metadata=adata.obs[[sample_key,group_key]].drop_duplicates()
    if metadata[sample_key].duplicated().any():raise ValueError('A sample cannot belong to multiple groups.')
    metadata=metadata.set_index(sample_key).loc[samples].rename(columns={group_key:'group_id'})
    counts={};n_cells=pd.DataFrame(0,index=clusters,columns=samples,dtype=int)
    for cluster in clusters:
        mask=(adata.obs[cluster_key]==cluster).to_numpy()
        cells=np.flatnonzero(mask)
        codes=pd.Index(samples).get_indexer(adata.obs.loc[mask,sample_key])
        membership=sparse.csr_matrix((np.ones(len(cells)),(codes,np.arange(len(cells)))),shape=(len(samples),len(cells)))
        result=membership@x[cells]
        result=result.toarray() if sparse.issparse(result) else np.asarray(result)
        counts[cluster]=pd.DataFrame(result.T,index=adata.var_names,columns=samples)
        n_cells.loc[cluster]=np.bincount(codes,minlength=len(samples))
    return Pseudobulk(counts,n_cells,metadata)


def _bh(p):
    p=np.asarray(p);n=len(p);order=np.argsort(p,kind='stable');out=np.empty(n)
    out[order]=np.minimum(np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1],1)
    return out


def pb_ds(pb, design, contrasts, *, method='edgeR', min_cells=10, filter='both'):
    """Pseudobulk differential state tests, with local and per-contrast global BH.

    design is samples x coefficients; contrasts is coefficients x comparisons.
    Both must be named DataFrames. Requires >=2 biological samples per group
    after filtering and residual degrees of freedom. Unsupported cell-level
    mixed models, treat tests and limma-trend are deliberately excluded.
    Returns result tables and explicit reasons for omitted clusters.
    """
    if method not in ('edgeR','limma-voom') or filter not in ('both','genes','samples','none'):
        raise ValueError('Unsupported method/filter.')
    if not isinstance(min_cells,int) or isinstance(min_cells,bool) or min_cells<1:
        raise ValueError('min_cells must be a positive integer.')
    if (not isinstance(design,pd.DataFrame) or not design.index.equals(pb.samples.index)
            or not design.columns.is_unique or not np.isfinite(design.to_numpy()).all()):
        raise ValueError('design must be finite and indexed in pseudobulk sample order.')
    if (not isinstance(contrasts,pd.DataFrame) or not contrasts.index.equals(design.columns)
            or not contrasts.columns.is_unique or not len(contrasts.columns)
            or not np.isfinite(contrasts.to_numpy()).all() or np.any((contrasts==0).all(axis=0))):
        raise ValueError('contrasts must align with design coefficients and be nonzero.')
    import pyedger as edge
    tables=[];skipped={}
    groups=pd.unique(pb.samples.group_id)
    for cluster,original in pb.counts.items():
        keep=pb.n_cells.loc[cluster]>=min_cells
        y=original.loc[:,keep];d=design.loc[y.columns]
        if filter in ('samples','both') and y.shape[1]:
            libs=y.sum(axis=0).to_numpy()
            log=np.log(libs,where=libs>0,out=np.full_like(libs,-np.inf))
            finite=log[np.isfinite(log)]
            cutoff=np.median(finite)-3*1.4826*np.median(np.abs(finite-np.median(finite))) if len(finite) else np.inf
            y=y.loc[:,log>=cutoff];d=d.loc[y.columns]
        counts=pb.samples.loc[y.columns,'group_id'].value_counts().reindex(groups,fill_value=0)
        if (counts<2).any() or len(d)<=d.shape[1] or np.linalg.matrix_rank(d)<d.shape[1]:
            skipped[cluster]='Insufficient group replication or rank-deficient/saturated design.';continue
        y=y.loc[y.sum(axis=1)!=0]
        if filter in ('genes','both') and y.size and y.to_numpy().max()>100:
            y=y.loc[np.asarray(edge.filterByExpr(y,design=d),dtype=bool)]
        y=y.loc[:,(y>0).any(axis=0)];d=d.loc[y.columns]
        remaining=pb.samples.loc[y.columns,'group_id'].value_counts().reindex(groups,fill_value=0)
        if not len(y) or (remaining<2).any() or len(d)<=d.shape[1] or np.linalg.matrix_rank(d)<d.shape[1]:
            skipped[cluster]='No estimable model after gene/empty-library filtering.';continue
        obj=edge.calcNormFactors(edge.DGEList(y))
        if method=='edgeR':
            obj=edge.estimateDisp(obj,design=d)
            # edgeR 4.8 glmQLFit.DGEList uses the mean trended dispersion of
            # the 10% most abundant genes when estimateDisp has already run.
            # pyedger 0.2.3 otherwise enters the no-dispersion Cox-Reid path.
            top=np.argsort(-np.asarray(obj.AveLogCPM),kind='stable')[:int(np.ceil(.1*len(y))) ]
            dispersion=float(np.mean(np.asarray(obj.trended_dispersion)[top]))
            fit=edge.glmQLFit(obj,design=d,dispersion=dispersion)
        else:
            import pylimma as limma
            lib=obj.samples['lib.size'].to_numpy()*obj.samples['norm.factors'].to_numpy()
            v=limma.voom(y,design=d,lib_size=lib)
            # Preserve muscat's explicit n_cells weights argument to lmFit.
            fit=limma.lmFit(v,design=d,weights=pb.n_cells.loc[cluster,y.columns].to_numpy())
        for comparison in contrasts.columns:
            if method=='edgeR':
                result=edge.glmQLFTest(fit,contrast=contrasts[comparison].to_numpy())
                table=edge.topTags(result,n=np.inf,sort_by='none').copy()
                table=table.rename(columns={'PValue':'p_val','FDR':'p_adj.loc'})
            else:
                model=limma.contrasts_fit(fit,contrasts[[comparison]])
                model=limma.eBayes(model,trend=False,robust=False)
                table=limma.topTable(model,number=np.inf,sort_by='none').copy()
                table=table.rename(columns={'P.Value':'p_val','adj.P.Val':'p_adj.loc'})
            table['gene']=y.index.to_numpy();table['cluster_id']=cluster;table['contrast']=comparison
            tables.append(table.reset_index(drop=True))
    if not tables:raise ValueError(f'No clusters are testable: {skipped}')
    table=pd.concat(tables,ignore_index=True)
    for comparison,indices in table.groupby('contrast',sort=False).groups.items():
        table.loc[indices,'p_adj.glb']=_bh(table.loc[indices,'p_val'])
    return {'table':table,'skipped':skipped,'design':design.copy(),'contrasts':contrasts.copy()}
