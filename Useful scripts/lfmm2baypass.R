lfmm2baypass<-function(geno.lfmm.file,metadata.lfmm.file,minpopsize=0,outprefix=""){
 geno.lfmm=read.table(geno.lfmm.file)
 metadata.lfmm=read.table(metadata.lfmm.file,h=T,row.names=1)
 outgenofile="geno.baypass" ; outcovfile="popcov.baypass"      
  if(nchar(outprefix)>0){
    outgenofile=paste(outprefix,outgenofile,sep=".")
    outcovfile=paste(outprefix,outcovfile,sep=".")
  }
 
 YYalt=t(rowsum(geno.lfmm,metadata.lfmm$pop))
 YYref=t(rowsum(geno.lfmm*0+2,metadata.lfmm$pop)) - YYalt
 tmp.npop=ncol(YYalt) ; tmp.nsnp=nrow(YYalt)
 genobaypass=matrix(0,tmp.nsnp,tmp.npop*2)
 tmp.idx=(1:tmp.npop)*2
 genobaypass[,tmp.idx-1]=YYalt
 rm(YYalt)
 genobaypass[,tmp.idx]=YYref
 rm(YYref)
 env.mean=t(rowsum(metadata.lfmm[,-1],metadata.lfmm$pop)/rowsum(metadata.lfmm[,-1]*0+1,metadata.lfmm$pop))
 
 if(minpopsize>0){
   popsize=table(metadata.lfmm$pop)
   pop.sel=which(popsize>=minpopsize)
   env.mean=env.mean[,pop.sel]
   genobaypass=genobaypass[,sort(c(pop.sel*2 - 1,pop.sel*2))]
 }
 
 write.table(file=outgenofile,genobaypass,quote=F,row.names=F,col.names=F)
 write.table(file=outcovfile,env.mean,quote=F,row.names=F,col.names=F)
}

#############################
#############################
# Exemples
#Full data
lfmm2baypass(geno.lfmm.file="~/Bureau/2022_SSMPG_BayPass/GIT_depot/SSMPG2022/Osuah tree/genotypes.lfmm",
             metadata.lfmm.file="~/Bureau/2022_SSMPG_BayPass/GIT_depot/SSMPG2022/Osuah tree/metadata.txt")
#selecting only a subset of the largest pop (with at least 6 sampled individual)
lfmm2baypass(geno.lfmm.file="~/Bureau/2022_SSMPG_BayPass/GIT_depot/SSMPG2022/Osuah tree/genotypes.lfmm",
             metadata.lfmm.file="~/Bureau/2022_SSMPG_BayPass/GIT_depot/SSMPG2022/Osuah tree/metadata.txt",
             minpopsize = 6,outprefix = "minpop6")

