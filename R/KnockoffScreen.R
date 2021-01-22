
KS.prelim<-function(Y, X=NULL, id=NULL, out_type="C", B=1000){
  ##Preliminary
  Y<-as.matrix(Y);n<-nrow(Y)

  if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
  X0<-cbind(rep(1,n),X0)

  if(out_type=="C"){nullglm<-glm(Y~0+X0,family=gaussian)}
  if(out_type=="D"){nullglm<-glm(Y~0+X0,family=binomial)}

  if (length(id)==0){id<-1:n}

  mu<-nullglm$fitted.values;Y.res<-Y-mu
  #permute the residuals
  index<-sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res<-Y.res[as.vector(index)]
  re.Y.res<-matrix(temp.Y.res,length(Y),B)

  #prepare invserse matrix for covariates
  if(out_type=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),length(Y))}
  inv.X0<-solve(t(X0)%*%(v*X0))

  #prepare the preliminary features
  result.prelim<-list(Y=Y,id=id,n=n,X0=X0,nullglm=nullglm,out_type=out_type,re.Y.res=re.Y.res,inv.X0=inv.X0)
  return(result.prelim)
}

KS.VCF.chr<-function(result.prelim,seq.filename,window.bed,M=5,thres.single=0.01,CADD.filename=NULL,GenoNet.filename=NULL,midout.dir=NULL,temp.dir=NULL,jobtitle=NULL,Gsub.id=NULL,impute.method='fixed'){

  time.start<-proc.time()[3]
  region.step=0.1*10^6

  old <- options()         # code line i
  on.exit(options(old)) 	 # code line i+1
  options(scipen = 999) #remove scientific notations
  chr<-window.bed[1,1]
  pos.min<-min(as.numeric(window.bed[,2]))
  pos.max<-max(as.numeric(window.bed[,3]))
  max.size<-max(as.numeric(window.bed[,3])-as.numeric(window.bed[,2]))

  start<-pos.min;null.model<-F;
  result.summary<-c();result.summary.single<-c()
  G.now<-c();G.before<-c();G.after<-c()
  while (start<pos.max){
    #if(start==200001){break}
    end<-min(pos.max,start+region.step)
    print(paste0(percent((start-pos.min)/(pos.max-pos.min)),' finished. Time used: ', round(proc.time()[3]-time.start,digits=1),'s. Scan ',start,'-',end))

    #8514
    #range<-paste0('chr',chr,":",max(start-0.1*10^6,1),"-",end+0.1*10^6)
    if(length(G.after)==0){
      range<-paste0(chr,":",max(start),"-",end)
      G.now <- t(readVCFToMatrixByRange(seq.filename, range,annoType='')[[1]])
      range<-paste0(chr,":",max(end,1),"-",end+region.step)
      G.after <- t(readVCFToMatrixByRange(seq.filename, range,annoType='')[[1]])
    }else{
      G.now <- G.after
      range<-paste0(chr,":",max(end,1),"-",end+region.step)
      G.after <- t(readVCFToMatrixByRange(seq.filename, range,annoType='')[[1]])
    }
    G<-cbind(G.before,G.now,G.after)
    G.before<-G.now

    #match phenotype id and genotype id
    if(length(Gsub.id)==0){match.index<-match(result.prelim$id,rownames(G))}else{
      match.index<-match(result.prelim$id,Gsub.id)
    }
    if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
    }
    G<-as.matrix(G[match.index,])
    if(length(G)==0){
      msg<-'Number of variants in the specified range is 0'
      warning(msg,call.=F)
      start<-end
      next
    }else{
      if(ncol(G)==1){
        msg<-'Number of variants in the specified range is 1'
        warning(msg,call.=F)
        start<-end
        next
      }}
    # missing genotype imputation
    G[G<0 | G==9]<-NA
    N_MISS<-sum(is.na(G))
    MISS.freq<-apply(is.na(G),2,mean)
    if(N_MISS>0){
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
      warning(msg,call.=F)
      G<-Impute(G,impute.method)
    }
    MAF<-apply(G,2,mean)/2
    G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
    MAF<-apply(G,2,mean)/2
    MAC<-apply(G,2,sum)
    s<-apply(G,2,sd)
    SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF) & MISS.freq<0.1)# & MAC>5

    pos<-as.numeric(gsub("^.*\\:","",colnames(G)))
    check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1 & pos>start & pos<end)
    if(length(check.index)<=1 ){
      msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
      warning(msg,call.=F)
      start<-end
      next
    }
    G<-Matrix(G[,SNP.index])
    MAF<-apply(G,2,mean)/2
    #get positions
    pos<-as.numeric(gsub("^.*\\:","",colnames(G)))
    #p.single<-Get.p(G,result.prelim)
    #G_k<-create.MK(G,pos,M=M,FA=-log10(p.single))
    #G_k<-create.MK(G,pos,M=M)
    G_k<-create.MK(G,pos,M=M,maxN.neighbor=Inf,maxBP.neighbor=0.1*10^6)

    ##single variant test
    single.index<-which(pos>start & pos<end & MAF>thres.single)
    if(length(single.index)>=1){
      G.single<-matrix(NA,nrow(G),length(single.index))
      G.single_k<-array(NA,dim=c(dim(G_k)[1],dim(G_k)[2],length(single.index)))
      G.single<-as.matrix(G[,single.index]);colnames(G.single)<-colnames(G)[single.index]
      G.single_k[,,]<-G_k[,,single.index]
      p.single.single<-Get.p(G.single,result.prelim)
      p.single_k.single<-apply(G.single_k,1,function(x) Get.p(x,result.prelim=result.prelim))
      if(ncol(G.single)==1){p.single_k.single<-matrix(p.single_k.single,1,dim(G_k)[1])}

      W<-(-log10(p.single.single)-apply(-log10(p.single_k.single),1,median))*(-log10(p.single.single)>=apply(-log10(p.single_k.single),1,max))
      W[is.na(W)]<-0
      MK.stat<-MK.statistic(-log10(p.single.single),-log10(p.single_k.single),method='median')
      temp.summary.single<-cbind(chr,pos[single.index],pos[single.index],
                                 pos[single.index],pos[single.index],
                                 MK.stat,
                                 W,p.single.single,
                                 p.single_k.single,MAF[single.index])
      colnames(temp.summary.single)<-c('chr','start','end','actual_start','actual_end',
                                       'kappa','tau',
                                       'W_KS','P_KS',paste0('P_KS_k',1:M),'MAF')
      if(length(midout.dir)!=0){
        write.table(temp.summary.single,paste0(midout.dir,jobtitle,'_single_',chr,':',start,'-',end,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
      }
      result.summary.single<-rbind(result.summary.single,temp.summary.single)
    }

    ##rare variants
    index<-which(pos>start & pos<end)
    G.window<-matrix(NA,nrow(G),length(index))
    G.window_k<-array(NA,dim=c(dim(G_k)[1],dim(G_k)[2],length(index)))
    G.window<-as.matrix(G[,index]);colnames(G.window)<-colnames(G)[index]
    G.window_k[,,]<-G_k[,,index]

    MAF<-apply(G.window,2,mean)/2
    MAC<-apply(G.window,2,sum)
    pos<-as.numeric(gsub("^.*\\:","",colnames(G.window)))

    window.temp<-window.bed[window.bed[,2]<max(pos) & window.bed[,3]>min(pos),]
    if(length(nrow(window.temp))==0){window.temp<-matrix(window.temp,1,ncol(window.bed))}
    window.matrix0<-matrix(apply(window.temp,1,function(x)as.numeric(pos>=x[2] & pos<x[3])),length(pos),nrow(window.temp))
    window.string<-apply(window.matrix0,2,function(x)paste(as.character(x),collapse = ""))

    window.MAC<-apply(MAC*window.matrix0,2,sum)
    window.index<-intersect(match(unique(window.string),window.string),which(apply(window.matrix0,2,sum)>1 & window.MAC>=10))
    if(length(window.index)==0){
      start<-end
      next
    }
    window.matrix0<-as.matrix(window.matrix0[,window.index])

    window.matrix<-Matrix(window.matrix0)
    window.summary<-cbind(window.temp[window.index,2],window.temp[window.index,3],t(apply(window.matrix,2,function(x)c(min(pos[which(x==1)]),max(pos[which(x==1)])))))

    weight.beta<-dbeta(apply(G.window,2,mean)/2,1,25)
    weight.matrix<-cbind(MAC<=10,(MAF<0.01&MAC>5)*weight.beta,(MAF>=0.01)*weight.beta)
    colnames(weight.matrix)<-c('MAC<=10','MAF<0.01&MAC>5&Beta','MAF>=0.01Beta')
    #CADD.filename<-'/oak/stanford/groups/zihuai/GeneticsResources/CADD_hg38/whole_genome_SNVs.tsv.gz'
    if(length(CADD.filename)!=0){
      job_id<-paste0(chr,':',pos.min,'-',pos.max)
      dir.string<-paste0('tabix ',CADD.filename,' ',sub(".*chr", "", chr),':',min(pos),'-',max(pos),' > ',temp.dir,'Temp_FA_',jobtitle,'_',job_id,'.txt')
      system(dir.string)
      score<-try(data.frame(fread(paste0(temp.dir,'Temp_FA_',jobtitle,'_',job_id,'.txt'))),silent = T)
      #index<-match(pos,score[,2])-1
      #CADD<-(score[index+1,6]+score[index+2,6]+score[index+3,6])/3
      index<-rep(match(pos,score[,2])-1,each=3)+1:3
      temp.score<-score[index,]
      CADD<-as.matrix(tapply(temp.score[,6],temp.score[,2],mean))
      colnames(CADD)<-'MAF<0.01&MAC>5&CADD'
      weight.matrix<-cbind(weight.matrix,(MAF<0.01&MAC>5)*CADD)
      #weight.matrix<-Matrix(cbind(MAC<10,(MAF<0.01)*dbeta(apply(G,2,mean)/2,1,25),dbeta(apply(G,2,mean)/2,1,25),rep(1,length(MAF)),CADD))
      #colnames(weight.matrix)<-c('MAC<10','MAF<0.01&Beta','Beta','Equal','CADD')
    }
    #GenoNet.filename=paste0('/oak/stanford/groups/zihuai/GenoNet/GenoNetScores_byChr/GenoNet_',chr,'.bed.gz')
    if(length(GenoNet.filename)!=0){
      job_id<-paste0(chr,':',pos.min,'-',pos.max)
      dir.string<-paste0('tabix ',GenoNet.filename,' chr',sub(".*chr", "", chr),':',max(0,min(pos)-100),'-',max(pos)+100,' > ',temp.dir,'Temp_FA_',jobtitle,'_',job_id,'.txt')
      system(dir.string)
      score<-try(data.frame(fread(paste0(temp.dir,'Temp_FA_',jobtitle,'_',job_id,'.txt'))),silent = T)
      if(length(score)==0){score<-matrix(1,length(pos),131)}
      GenoNet.all<-score[sapply(pos,function(s)which.min(abs(s-score[,2]))),]
      if(length(nrow(GenoNet.all))==0){GenoNet.all<-matrix(GenoNet.all,1,131)}
      colnames(GenoNet.all)<-c('chr','start','end','name',
                               "E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009", "E010",
                               "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018", "E019", "E020",
                               "E021", "E022", "E023", "E024", "E025", "E026", "E027", "E028", "E029", "E030",
                               "E031", "E032", "E033", "E034", "E035", "E036", "E037", "E038", "E039", "E040",
                               "E041", "E042", "E043", "E044", "E045", "E046", "E047", "E048", "E049", "E050",
                               "E051", "E052", "E053", "E054", "E055", "E056", "E057", "E058", "E059", "E061",
                               "E062", "E063", "E065", "E066", "E067", "E068", "E069", "E070", "E071", "E072",
                               "E073", "E074", "E075", "E076", "E077", "E078", "E079", "E080", "E081", "E082",
                               "E083", "E084", "E085", "E086", "E087", "E088", "E089", "E090", "E091", "E092",
                               "E093", "E094", "E095", "E096", "E097", "E098", "E099", "E100", "E101", "E102",
                               "E103", "E104", "E105", "E106", "E107", "E108", "E109", "E110", "E111", "E112",
                               "E113", "E114", "E115", "E116", "E117", "E118", "E119", "E120", "E121", "E122",
                               "E123", "E124", "E125", "E126", "E127", "E128", "E129")

      StemCell<-c('E022','E007','E004','E002','E021','E009','E010','E001','E015','E018','E016','E020','E014','E019','E024','E008','E003','E012','E011')
      Blood<-c('E115', 'E123','E030','E029','E124','E035','E036','E051','E050','E034','E046','E041','E047','E048','E038','E045','E044','E043','E039','E042','E040','E037','E112','E093','E062','E033','E116','E031','E032')
      ConnectiveTissue<-c('E122','E120','E121','E025','E023','E049','E026','E129','E126','E052','E125','E055','E056','E017','E128','E114','E117','E028','E057','E058','E119','E127')
      Brain<-c('E071','E074','E073','E068','E067','E069','E071','E072')
      InternalOrgans<-c('E027','E059','E061','E065','E097','E086','E087','E100','E105','E104','E095','E096','E113','E079','E094','E098')
      FetalBrain<-c('E081','E070','E082','E054','E053')
      FetalTissue1<-c('E005','E099','E013','E006','E083')
      Muscle<-c('E108','E107','E063','E078','E103','E076','E111')
      FetalTissue2<-c('E091','E092','E089','E090','E088','E080')
      GI<-c('E066','E110','E109','E106','E075','E077','E101','E102','E118','E085','E084')

      if(nrow(GenoNet.all)==1){
        GenoNet<-cbind(max(GenoNet.all[,StemCell]),
                       max(GenoNet.all[,Blood]),
                       max(GenoNet.all[,ConnectiveTissue]),
                       max(GenoNet.all[,Brain]),
                       max(GenoNet.all[,InternalOrgans]),
                       max(GenoNet.all[,FetalBrain]),
                       max(GenoNet.all[,FetalTissue1]),
                       max(GenoNet.all[,Muscle]),
                       max(GenoNet.all[,FetalTissue2]),
                       max(GenoNet.all[,GI]))
      }else{
        GenoNet<-cbind(apply(GenoNet.all[,StemCell],1,max),
                       apply(GenoNet.all[,Blood],1,max),
                       apply(GenoNet.all[,ConnectiveTissue],1,max),
                       apply(GenoNet.all[,Brain],1,max),
                       apply(GenoNet.all[,InternalOrgans],1,max),
                       apply(GenoNet.all[,FetalBrain],1,max),
                       apply(GenoNet.all[,FetalTissue1],1,max),
                       apply(GenoNet.all[,Muscle],1,max),
                       apply(GenoNet.all[,FetalTissue2],1,max),
                       apply(GenoNet.all[,GI],1,max))
      }

      colnames(GenoNet)<-paste0('MAF<0.01&MAC>5&',c('StemCell','Blood','ConnectiveTissue','Brain','InternalOrgans','FetalBrain','FetalTissue1','Muscle','FetalTissue2','Gastrointestinal'))
      weight.matrix<-cbind(weight.matrix,(MAF<0.01&MAC>5)*GenoNet)
    }
    weight.matrix<-Matrix(weight.matrix)

    KS.fit<-KS.test(G.window,G.window_k,result.prelim,window.matrix=window.matrix,weight.matrix=weight.matrix)
    p.A<-KS.fit$p.KS;p.A_k<-KS.fit$p.KS_k
    p.A.common<-KS.fit$p.KS.common;p.A_k.common<-KS.fit$p.KS_k.common
    p.A.rare<-KS.fit$p.KS.rare;p.A_k.rare<-KS.fit$p.KS_k.rare

    p.individual<-KS.fit$p.individual
    #Knockoff statistics
    W<-(-log10(p.A)-apply(-log10(p.A_k),1,median))*(-log10(p.A)>=apply(-log10(p.A_k),1,max))
    W.common<-(-log10(p.A.common)-apply(-log10(p.A_k.common),1,median))*(-log10(p.A.common)>=apply(-log10(p.A_k.common),1,max))
    W.rare<-(-log10(p.A.rare)-apply(-log10(p.A_k.rare),1,median))*(-log10(p.A.rare)>=apply(-log10(p.A_k.rare),1,max))
    #W[is.na(W)]<-W.common[is.na(W.common)]<-W.rare[is.na(W.rare)]<-0
    MK.stat<-MK.statistic(-log10(p.A),-log10(p.A_k),method='median')

    temp.summary<-cbind(chr,window.summary,
                        MK.stat,
                        W,p.A,p.A_k,
                        W.common,W.rare,
                        p.A.common,p.A.rare,
                        p.A_k.common,p.A_k.rare,
                        p.individual)
    colnames(temp.summary)[1:(7+3*(2+M))]<-c('chr','start','end','actual_start','actual_end',
                                             'kappa','tau',
                                             'W_KS','P_KS',paste0('P_KS_k',1:M),
                                             'W_KS_common','W_KS_rare',
                                             'P_KS_common','P_KS_rare',
                                             paste0('P_KS_k_common',1:M),paste0('P_KS_k_rare',1:M))

    if(length(midout.dir)!=0){
      write.table(temp.summary,paste0(midout.dir,jobtitle,'_window_',chr,':',start,'-',end,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
    }

    result.summary<-rbind(result.summary,temp.summary)

    start<-end
  }

  return(list(result.window=result.summary,result.single=result.summary.single))
}



KS.summary<-function(result.window,result.single,M){

  temp<-matrix(NA,nrow(result.single),ncol(result.window))
  colnames(temp)<-colnames(result.window)
  temp[,colnames(result.single)[-ncol(result.single)]]<-result.single[,colnames(result.single)[-ncol(result.single)]]

  result<-rbind(result.window,temp)
  result<-result[order(result[,2]),]
  result<-result[order(result[,1]),]

  q<-MK.q.byStat(result[,'kappa'],result[,'tau'],M=5)
  result.summary<-cbind(result[,1:5],q,result[,-(1:5)])
  colnames(result.summary)[6]<-'Qvalue'

  return(result.summary[,1:grep('P_KS_rare',colnames(result.summary))])
}


KS.test<-function(G,G_k,result.prelim,window.matrix,weight.matrix){
  G<-Matrix(G)
  MAF<-apply(G,2,mean)/2;MAC<-apply(G,2,sum)
  window.matrix<-Matrix(window.matrix);weight.matrix<-Matrix(weight.matrix)

  mu<-result.prelim$nullglm$fitted.values;
  Y.res<-result.prelim$Y-mu;re.Y.res<-result.prelim$re.Y.res
  X0<-result.prelim$X0;outcome<-result.prelim$out_type

  #Burden test
  p.burden<-matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  p.burden_k<-array(NA,dim=c(dim(G_k)[1],ncol(window.matrix),ncol(weight.matrix)))
  for (k in 1:ncol(weight.matrix)){
    temp.window.matrix<-weight.matrix[,k]*window.matrix
    X<-as.matrix(G%*%temp.window.matrix)
    X_k<-apply(G_k,1,function(x) x%*%temp.window.matrix)
    p.burden[,k]<-Get.p.base(X,result.prelim)
    p.burden_k[,,k]<-t(sapply(X_k,Get.p.base,result.prelim=result.prelim))
  }

  #Dispersion test
  p.dispersion<-matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  p.dispersion_k<-array(NA,dim=c(dim(G_k)[1],ncol(window.matrix),ncol(weight.matrix)))

  #proc.time()
  if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(G))}
  A<-t(G)%*%(v*G)
  B<-t(G)%*%(v*X0)
  C<-solve(t(X0)%*%(v*X0))
  K<-A-B%*%C%*%t(B)

  score<-t(G)%*%Y.res;re.score<-t(t(G)%*%re.Y.res)
  score_k<-sapply(1:dim(G_k)[1],function(s){t(G_k[s,,])%*%Y.res});
  re.score_k<-lapply(1:dim(G_k)[1],function(s) t(t(G_k[s,,])%*%re.Y.res))

  for (k in 1:ncol(weight.matrix)){
    #print(k)
    p.dispersion[,k]<-Get.p.SKAT(score,re.score,K,window.matrix,weight=(MAC>=10)*weight.matrix[,k],result.prelim)
    p.dispersion_k[,,k]<-t(sapply(1:dim(G_k)[1],function(s){Get.p.SKAT(score_k[,s],re.score_k[[s]],K,window.matrix,weight=(MAC>=10)*weight.matrix[,k],result.prelim)}))
  }
  #proc.time()

  p.single<-Get.p(G,result.prelim)
  p.single_k<-apply(G_k,1,function(x) Get.p(x,result.prelim=result.prelim))
  if(ncol(G)==1){p.single_k<-matrix(p.single_k,1,dim(G_k)[1])}
  p.V1<-Get.cauchy.scan(p.single,(MAC>=10 & MAF<0.01)*window.matrix)
  p.V1_k<-apply(p.single_k,2,Get.cauchy.scan,window.matrix=(MAC>=10 & MAF<0.01)*window.matrix)
  p.V2<-Get.cauchy.scan(p.single,(MAF>=0.01)*window.matrix)
  p.V2_k<-apply(p.single_k,2,Get.cauchy.scan,window.matrix=(MAF>=0.01)*window.matrix)
  if(ncol(window.matrix)==1){p.V1_k<-matrix(p.V1_k,1,dim(G_k)[1]);p.V2_k<-matrix(p.V2_k,1,dim(G_k)[1])}

  p.individual<-cbind(p.burden,p.dispersion,p.V1,p.V2);
  colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<0.01&MAC>5','singleCauchy_MAF>=0.01')

  p.KS<-as.matrix(apply(p.individual,1,Get.cauchy))
  p.KS_k<-matrix(sapply(1:dim(G_k)[1],function(s){apply(cbind(matrix(p.burden_k[s,,],dim(p.burden_k)[2],dim(p.burden_k)[3]),
                                                                matrix(p.dispersion_k[s,,],dim(p.dispersion_k)[2],dim(p.dispersion_k)[3]),
                                                                p.V1_k[,s],p.V2_k[,s]),1,Get.cauchy)}),dim(p.burden_k)[2],dim(p.burden_k)[1])
  test.common<-grep('MAF>=0.01',colnames(p.individual))
  p.KS.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
  p.KS_k.common<-matrix(sapply(1:dim(G_k)[1],function(s){apply(cbind(matrix(p.burden_k[s,,],dim(p.burden_k)[2],dim(p.burden_k)[3]),
                                                              matrix(p.dispersion_k[s,,],dim(p.dispersion_k)[2],dim(p.dispersion_k)[3]),
                                                              p.V1_k[,s],p.V2_k[,s])[,test.common,drop=FALSE],1,Get.cauchy)}),dim(p.burden_k)[2],dim(p.burden_k)[1])

  p.KS.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
  p.KS_k.rare<-matrix(sapply(1:dim(G_k)[1],function(s){apply(cbind(matrix(p.burden_k[s,,],dim(p.burden_k)[2],dim(p.burden_k)[3]),
                                                                     matrix(p.dispersion_k[s,,],dim(p.dispersion_k)[2],dim(p.dispersion_k)[3]),
                                                                     p.V1_k[,s],p.V2_k[,s])[,-test.common,drop=FALSE],1,Get.cauchy)}),dim(p.burden_k)[2],dim(p.burden_k)[1])

  return(list(p.KS=p.KS,p.KS_k=p.KS_k,p.KS.common=p.KS.common,p.KS_k.common=p.KS_k.common,p.KS.rare=p.KS.rare,p.KS_k.rare=p.KS_k.rare,p.individual=p.individual,p.single=p.single,p.single_k=p.single_k))
}

Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get.p.SKAT<-function(score,re.score,K,window.matrix,weight,result.prelim){

  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  X0<-result.prelim$X0;outcome<-result.prelim$out_type

  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2)
  K.temp<-weight*t(weight*K)

  #fast implementation by resampling based moment matching
  p0<-Get.p.moment(as.vector(t(score^2)%*%(weight*window.matrix)^2),re.score^2%*%(weight*window.matrix)^2)
  p<-p0#c()
  for(i in which(p0<0.01 |p0>=1)){
    #print(i)
    temp<-K.temp[window.matrix[,i]!=0,window.matrix[,i]!=0]
    if(sum(temp^2)==0){p[i]<-NA;next}

    #proc.time()
    lambda=eigen(temp,symmetric=T,only.values=T)$values
    temp.p<-davies(Q[i],lambda,acc=10^(-6))$Qq
    #temp.p<-SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
    #proc.time()

    #proc.time()
    #temp.p<-Get_Liu_PVal.MOD(Q[i],temp*2)$p.value
    #proc.time()

    #temp.p<-davies(Q[i],lambda,acc=10^(-6))$Qq
    if(temp.p > 1 || temp.p <= 0 ){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
    }
    p[i]<-temp.p
  }
  return(as.matrix(p))
}


#percentage notation
percent <- function(x, digits = 3, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}



create.MK <- function(X,pos,M=5,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=100000) {

  sparse.fit<-sparse.cor(X)
  cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  #temp.M<-matrix(1,M+1,M+1)
  #cov.M<-kronecker(temp.M,cov.X)

  #X_k<-matrix(0,nrow(X),ncol(X));index.exist<-c()
  X_k<-array(0,dim=c(M,nrow(X),ncol(X)));index.exist<-c()
  for (k in unique(clusters)){
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)){
      #print(i)
      index.pos<-which(pos>=max(pos[i]-maxBP.neighbor,pos[1]) & pos<=min(pos[i]+maxBP.neighbor,pos[length(pos)]))
      temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0

      index<-order(temp,decreasing=T)
      index<-setdiff(index[1:min(length(index),sum(temp>0.05),floor((nrow(X))^(1/3)),maxN.neighbor)],i)

      y<-X[,i]
      if(length(index)==0){fitted.values<-mean(y)}else{

        x<-X[,index,drop=F];temp.xy<-rbind(mean(y),crossprod(x,y)/length(y)-colMeans(x)*mean(y))
        x.exist<-c()
        for(j in 1:M){
          x.exist<-cbind(x.exist,X_k[j,,intersect(index,index.exist)])
        }
        temp.xy<-rbind(temp.xy,crossprod(x.exist,y)/length(y)-colMeans(x.exist)*mean(y))

        temp.cov.cross<-sparse.cov.cross(x,x.exist)$cov
        temp.cov<-sparse.cor(x.exist)$cov
        temp.xx<-cov.X[index,index]
        temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))

        temp.xx<-cbind(0,temp.xx)
        temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

        pca.fit<-princomp(covmat=temp.xx)
        v<-pca.fit$loadings
        cump<-cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc<-which(cump>=0.999)[1]#nrow(temp.xx)#nrow(temp.xx)#
        pca.index<-intersect(1:n.pc,which(pca.fit$sdev!=0))#which(cump<=0.99)
        #calculate
        #inverse ZZ matrix
        temp.inv<-v[,pca.index,drop=F]%*%(pca.fit$sdev[pca.index]^(-2)*t(v[,pca.index,drop=F]))
        #beta coefficients
        temp.beta<-temp.inv%*%temp.xy

        temp.j<-1
        fitted.values<-temp.beta[1]+crossprod(t(x),temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
        temp.j<-temp.j+ncol(x)
        for(j in 1:M){
          temp.x<-as.matrix(X_k[j,,intersect(index,index.exist)])
          if(ncol(temp.x)>=1){
            fitted.values<-fitted.values+crossprod(t(temp.x),temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
          }
          temp.j<-temp.j+ncol(temp.x)
        }
      }
      residuals<-y-fitted.values
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)

      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[j,,which(clusters==k)]<-cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F]
    }
  }
  return(X_k)
}

max.nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}

Get.p.base<-function(X,result.prelim){
  X<-as.matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  outcome<-result.prelim$out_type
  if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(X))}
  p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
  p[is.na(p)]<-NA
  return(p)
}

Get.p<-function(X,result.prelim){
  X<-as.matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  outcome<-result.prelim$out_type
  if(outcome=='D'){
    p<-ScoreTest_SPA(t(X),result.prelim$Y,result.prelim$X,method=c("fastSPA"),minmac=-Inf)$p.value
  }else{
    v<-rep(as.numeric(var(Y.res)),nrow(X))
    p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
  }
  return(as.matrix(p))
}

MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  kappa<-apply(T.temp,1,which.max)-1
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}

MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
    return(tau[b][ok[length(ok)]])
  }else{return(Inf)}
}

MK.threshold<-function (T_0,T_k, fdr = 0.1,method='median',Rej.Bound=10000){
  stat<-MK.statistic(T_0,T_k,method=method)
  kappa<-stat[,1];tau<-stat[,2]
  t<-MK.threshold.byStat(kappa,tau,M=ncol(T_k),fdr=fdr,Rej.Bound=Rej.Bound)
  return(t)
}


MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau))
  for(i in 1:length(b)){
    q[b[i]]<-min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}


Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))
  #window.matrix.MAC10<-(MAC>=10)*window.matrix0

  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  #cct.stat<-as.numeric(t(temp)%*%window.matrix.MAC10/apply(window.matrix.MAC10,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
}

Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean<-apply(re.Q,2,mean)
  re.variance<-apply(re.Q,2,var)
  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
  #re.p[re.p==1]<-0.99
  return(re.p)
}


Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}


