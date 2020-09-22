X<-read.table(file='KDH/KDH13toEL10_75_besthits_order.txt')
XX<-read.table('KDH/gene_pos1.txt')

i<-unique(X[,13])[1]
Y<-X[X[,13]==i,]
I<-unique(Y[,1])
I
I1<-unique(Y[,16])


Final1<-c()
for (i in unique(X[,13])){
  PepDat<-c()
  Y<-X[X[,13]==i,]
  I<-unique(Y[,1])
  Mat<-c()
  for (I2 in I){
    line<-Y[Y[,1]==I2,]
    INDEX<-sort(line[,11],decreasing = F,index.return = T)
    line<-line[(INDEX$ix)[1],]
    PepDat<-rbind(PepDat,line[1,])
  }
  ELChr<-unique(PepDat[,16][1])
  slope<-PepDat[PepDat[,16]==ELChr,][,17]
  PepDat[PepDat[,16]==ELChr,][,2]
  for (ge in unique(gsub(pattern = '.1',replacement = '',x=PepDat[PepDat[,16]==ELChr,][,2],fixed=TRUE))){
    Maxi<-c()
    Mini<-c()
    Maxi<-c(Maxi,XX[XX[,4]==ge,][,3])
    Mini<-c(Mini,XX[XX[,4]==ge,][,2])
    CHRO<-unique(XX[XX[,4]==ge,1])
    if (length(Mini) == 0){
      Mini<-c('1')
      Maxi<-c('1')
      CHRO<-'NA'
    }
  }
  Maxmax<-max(Maxi)
  Minmin<-min(Mini)
  Means<-mean(PepDat[PepDat[,16]==ELChr,][,11])
  Mat<-cbind(seq(1,length(slope)),slope)
  Mat<-as.matrix(Mat)
  #plot(Mat)
  LM<-lm.fit(Mat,Mat)
  Final<-cbind(as.character(ELChr),LM$coefficients[3],as.character(i),median(slope),length(slope),Means,
               Maxmax,Minmin,CHRO)
  Final1<-rbind(Final1,Final)
  #print(LM$coefficients[3])
  #print(i)
  #print(median(slope))
}

Final1<-as.matrix(Final1)
colnames(Final1)<-c('chr','slope','contig','pos','len','means','max','min','CHR')

Final2<-c()
for (i in sort(unique(Final1[,1]))){
  SortedFinal<-Final1[Final1[,1]==i,]
  INDEX1<-sort(SortedFinal[,4],decreasing = F,index.return = T)
  print(SortedFinal[INDEX1$ix,])
  Final2<-rbind(Final2,SortedFinal[INDEX1$ix,])
}

Final2<-rbind(Final2,Final1[Final1[,1]=='EL10As22',])
Final2
write.table(Final2,file = 'KDH/Gene_order.txt')

ChrOrderMatrix<-read.table('KDH/Gene_order.txt')
ChrOrderMatrix$CHR<-paste('Chr',ChrOrderMatrix$CHR,sep = '')
write.table(ChrOrderMatrix,file = 'KDH/Gene_order1.txt')
