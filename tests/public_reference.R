# Public Kang et al. IFN-stimulation subset distributed by muscat.
# Run after tests/reference.R with the same pinned source checkout.
source('tests/reference.R')
load(file.path(ref,'data','example_sce.rda'))
sce<-example_sce[seq_len(200),]
pb<-aggregateData(sce,verbose=FALSE)
d<-model.matrix(~group_id,data=as.data.frame(colData(pb)))
rownames(d)<-colnames(pb)
ct<-matrix(c(0,1),2,1,dimnames=list(colnames(d),'stim'))
write.csv(as.matrix(counts(sce)),'data/public_counts.csv')
write.csv(as.data.frame(colData(sce)),'data/public_metadata.csv')
write.csv(d,'data/public_design.csv');write.csv(ct,'data/public_contrast.csv')
for (method in c('edgeR','limma-voom')) {
 result<-pbDS(pb,method=method,design=d,contrast=ct,min_cells=10,filter='both',verbose=FALSE)
 write.csv(do.call(rbind,result$table$stim),paste0('data/public_',method,'.csv'),row.names=FALSE)
}
