library(peer)
fn='RNAseq_healthy_controls.csv'
expr = read.csv(fn, header=T)

model = PEER()
PEER_setPhenoMean(model,as.matrix(expr[-1])) # Skip first column LABID
PEER_setNk(model,10) # Calcualte 10 PEER factors

PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model) # Get factors (NxK)
write.table(factors, file="PEER_control.txt")
