cuts <- function(x)
{
  n <- length(x) %/% 4
  map <- rep(c(rep(TRUE,4),FALSE), n)
  result <- rep(NA, n*5)
  result[map] <- x
  result
}

##Circos Plot
{
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  
  #install.packages('RCircos')
  library(RCircos)
  BeetCent<-read.table('/Volumes/UNTITLED/Projects/GSS/R_Circos/Tony_CHmarks/Paul_Blast/Centro.txt')
  
  BeetCent
  
  chr.exclude <- NULL;
  cyto.info <- BeetCent;
  tracks.inside <- 10;
  tracks.outside <- 0;
  RCircos.Set.Core.Components(cyto.info, chr.exclude,
                              tracks.inside, tracks.outside);
  
  RCircos.Set.Plot.Area()
  
  par(mai=c(0.25, 0.25, 0.25, 0.25));
  plot.new();
  plot.window(c(-2.5,2.5), c(-2.5, 2.5));
  #RCircos.Chromosome.Ideogram.Plot(tick.interval=2,ideo.width=3);
  #RCircos.Draw.Chromosome.Ideogram(ideo.pos=1.9, ideo.width=0.1)
  RCircos.Highligh.Chromosome.Ideogram(highlight.pos=2.3, highlight.width=0)
  RCircos.Ideogram.Tick.Plot(tick.interval=3, track.for.ticks=2.3)
  RCircos.Label.Chromosome.Names(chr.name.pos=NULL)
  
  source('/Volumes/UNTITLED/Projects/GSS/R_Circos/Tony_CHmarks/CYtoBand_fun.R')
  source('/Volumes/UNTITLED/Projects/GSS/R_Circos/Tony_CHmarks/CHRStain.R')
  
  GeneDen<-read.table('KDH/GeneDen.txt')
  RepeatDen<-read.table('KDH/RepeatDen.txt')
}

#GeneOrder - Un-Ordered
{
  PaulN<-1.85
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.07
  
  ChrOrderMatrix<-read.table('KDH/Gene_order.txt')
  ChrOrderMatrix$CHR<-paste('Chr',ChrOrderMatrix$CHR,sep = '')
  ChrOrderMatrix<-as.matrix(ChrOrderMatrix)
  orderChr1<-sort(ChrOrderMatrix[,3],decreasing = F,index.return = T)
  ChrOrderMatrix[orderChr1$ix,]
  CG<-c()
  for (i in ChrOrderMatrix[orderChr1$ix,9]){
    if (i == "Chr1"){
      CG<-c(CG,'orange')
    }
    if (i == "Chr2"){
      CG<-c(CG,'coral')
    }
    if (i == "Chr3"){
      CG<-c(CG,'brown')
    }
    if (i == "Chr4"){
      CG<-c(CG,'green3')
    }
    if (i=="Chr5"){
      CG<-c(CG,'coral4')
    }
    if (i == "Chr6"){
      CG<-c(CG,'red')
    }
    if (i == "Chr7"){
      CG<-c(CG,'darkgreen')
    }
    if (i == "Chr8"){
      CG<-c(CG,'magenta')
    }
    if (i == "Chr9"){
      CG<-c(CG,'pink1')
    }
  }
  
  for (P in seq(1,length(ChrOrderMatrix[,1]))){
    highlight.area=(ChrOrderMatrix[P,c(9,8,7)])
    #print(highlight.area)
    start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[2]));
    end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                  as.numeric(highlight.area[3]));
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    
    PolyX<-rbind(PolyX,c(polygon.x))
    PolyY<-rbind(PolyY,c(polygon.y))
    #print(N)
    N<-N+1
  }
  
  PolyX<-PolyX[-1,]
  PolyY<-PolyY[-1,]
  PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
  PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
  print('yes')
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=adjustcolor('white', alpha.f = 0.15))
  

  }

#GeneOrder - Ordered
{
  PaulN<-1.7
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.07
  
  ChrOrderMatrix<-read.table('KDH/Gene_order.txt')
  ChrOrderMatrix$CHR<-paste('Chr',ChrOrderMatrix$CHR,sep = '')
  ChrOrderMatrix<-as.matrix(ChrOrderMatrix)
  {
  CG<-c()
  for (i in ChrOrderMatrix[,9]){
    if (i == "Chr1"){
      CG<-c(CG,'orange')
    }
    if (i == "Chr2"){
      CG<-c(CG,'blue')
    }
    if (i == "Chr3"){
      CG<-c(CG,'brown')
    }
    if (i == "Chr4"){
      CG<-c(CG,'green3')
    }
    if (i=="Chr5"){
      CG<-c(CG,'coral4')
    }
    if (i == "Chr6"){
      CG<-c(CG,'red')
    }
    if (i == "Chr7"){
      CG<-c(CG,'darkgreen')
    }
    if (i == "Chr8"){
      CG<-c(CG,'magenta')
    }
    if (i == "Chr9"){
      CG<-c(CG,'pink1')
    }
  }
  }
  
  
  for (P in seq(1,length(ChrOrderMatrix[,1]))){
    highlight.area=(ChrOrderMatrix[P,c(9,8,7)])
    #print(highlight.area)
    start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[2]));
    end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                  as.numeric(highlight.area[3]));
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    
    PolyX<-rbind(PolyX,c(polygon.x))
    PolyY<-rbind(PolyY,c(polygon.y))
    #print(N)
    N<-N+1
  }
  
  PolyX<-PolyX[-1,]
  PolyY<-PolyY[-1,]
  PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
  PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
  print('yes')
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=adjustcolor('white', alpha.f = 0.1))


}  

#KDH13 - indels
{
PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH13UniqueIndels_for_Rcircos.txt")
PX<-c()
ssC<-PLCT[,1]
PLCT3<-c(0,0)
for (i in unique(PLCT[,1])){
  PLCT1<-PLCT[PLCT[,1]==i,3]
  HP<-hist(PLCT1,plot = F,breaks = 250 )
  length(HP$counts)
  length(HP$breaks)
  PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
  PLCT3<-rbind(PLCT3,PLCT2)
}
PLCT3<-PLCT3[-1,]

color.gradient <- function(x, colors=c("white","lightblue"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
CG<-color.gradient(HP$counts)
#plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')


setwd('/Volumes/UNTITLED/IMAD_2020/')
#POLYDRAW SNP

  N<-1
  PaulN<-1
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
    #POLYDRAW No DATA    
    {
      for (P in seq(1,length(PLCT3[,1]))){
        highlight.area=as.vector(PLCT3[P,c(1,2,2)])
        #print(highlight.area)
        start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                        as.numeric(highlight.area[2]));
        end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[3]));
        
        polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                       RCircos.Pos[end.pos:start.pos, 1]*innerPos);
        polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                       RCircos.Pos[end.pos:start.pos, 2]*innerPos);
        
        
        PolyX<-rbind(PolyX,c(polygon.x))
        PolyY<-rbind(PolyY,c(polygon.y))
        #print(N)
        N<-N+1
      }
      
      PolyX<-PolyX[-1,]
      PolyY<-PolyY[-1,]
      PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
      PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
      print('yes')
      polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
      PaulN<-PaulN+0.05
    }  

}
  
#KDH13 - SNP
{
  PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH13UniqueSNP_for_Rcircos.txt")
  PX<-c()
  ssC<-PLCT[,1]
  PLCT3<-c(0,0)
  for (i in unique(PLCT[,1])){
    PLCT1<-PLCT[PLCT[,1]==i,3]
    HP<-hist(PLCT1,plot = F,breaks = 250 )
    length(HP$counts)
    length(HP$breaks)
    PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
    PLCT3<-rbind(PLCT3,PLCT2)
  }
  PLCT3<-PLCT3[-1,]
  
  color.gradient <- function(x, colors=c("white","blue"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  CG<-color.gradient(HP$counts)
  #plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')
  
  
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  #POLYDRAW SNP
  
  N<-1
  PaulN<-1.06
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
  #POLYDRAW No DATA    
  {
    for (P in seq(1,length(PLCT3[,1]))){
      highlight.area=as.vector(PLCT3[P,c(1,2,2)])
      #print(highlight.area)
      start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[2]));
      end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[3]));
      
      polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 1]*innerPos);
      polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 2]*innerPos);
      
      
      PolyX<-rbind(PolyX,c(polygon.x))
      PolyY<-rbind(PolyY,c(polygon.y))
      #print(N)
      N<-N+1
    }
    
    PolyX<-PolyX[-1,]
    PolyY<-PolyY[-1,]
    PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
    PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
    print('yes')
    polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
    PaulN<-PaulN+0.05
  }  
  
}

#KDH19-17 - indel
{
  PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH19_17UniqueIndels_for_Rcircos.txt")
  PX<-c()
  ssC<-PLCT[,1]
  PLCT3<-c(0,0)
  for (i in unique(PLCT[,1])){
    PLCT1<-PLCT[PLCT[,1]==i,3]
    HP<-hist(PLCT1,plot = F,breaks = 250 )
    length(HP$counts)
    length(HP$breaks)
    PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
    PLCT3<-rbind(PLCT3,PLCT2)
  }
  PLCT3<-PLCT3[-1,]
  
  color.gradient <- function(x, colors=c("white","lightblue"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  CG<-color.gradient(HP$counts)
  #plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')
  
  
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  #POLYDRAW SNP
  
  N<-1
  PaulN<-1.2
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
  #POLYDRAW No DATA    
  {
    for (P in seq(1,length(PLCT3[,1]))){
      highlight.area=as.vector(PLCT3[P,c(1,2,2)])
      #print(highlight.area)
      start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[2]));
      end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[3]));
      
      polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 1]*innerPos);
      polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 2]*innerPos);
      
      
      PolyX<-rbind(PolyX,c(polygon.x))
      PolyY<-rbind(PolyY,c(polygon.y))
      #print(N)
      N<-N+1
    }
    
    PolyX<-PolyX[-1,]
    PolyY<-PolyY[-1,]
    PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
    PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
    print('yes')
    polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
    PaulN<-PaulN+0.05
  }  
  
}

#KDH19-17 - SNP
{
  PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH19_17UniqueSNP_for_Rcircos.txt")
  PX<-c()
  ssC<-PLCT[,1]
  PLCT3<-c(0,0)
  for (i in unique(PLCT[,1])){
    PLCT1<-PLCT[PLCT[,1]==i,3]
    HP<-hist(PLCT1,plot = F,breaks = 250 )
    length(HP$counts)
    length(HP$breaks)
    PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
    PLCT3<-rbind(PLCT3,PLCT2)
  }
  PLCT3<-PLCT3[-1,]
  
  color.gradient <- function(x, colors=c("white","red"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  CG<-color.gradient(HP$counts)
  #plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')
  
  
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  #POLYDRAW SNP
  
  N<-1
  PaulN<-1.25
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
  #POLYDRAW No DATA    
  {
    for (P in seq(1,length(PLCT3[,1]))){
      highlight.area=as.vector(PLCT3[P,c(1,2,2)])
      #print(highlight.area)
      start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[2]));
      end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[3]));
      
      polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 1]*innerPos);
      polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 2]*innerPos);
      
      
      PolyX<-rbind(PolyX,c(polygon.x))
      PolyY<-rbind(PolyY,c(polygon.y))
      #print(N)
      N<-N+1
    }
    
    PolyX<-PolyX[-1,]
    PolyY<-PolyY[-1,]
    PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
    PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
    print('yes')
    polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
    PaulN<-PaulN+0.05
  }  
  
}

#F1 - indel
{
  PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH_F1_P1_informative.vcf")
  PX<-c()
  ssC<-PLCT[,1]
  PLCT3<-c(0,0)
  for (i in unique(PLCT[,1])){
    PLCT1<-PLCT[PLCT[,1]==i,3]
    HP<-hist(PLCT1,plot = F,breaks = 250 )
    length(HP$counts)
    length(HP$breaks)
    PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
    PLCT3<-rbind(PLCT3,PLCT2)
  }
  PLCT3<-PLCT3[-1,]
  
  color.gradient <- function(x, colors=c("white","lightblue"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  CG<-color.gradient(HP$counts)
  #plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')
  
  
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  #POLYDRAW SNP
  
  N<-1
  PaulN<-1.3
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
  #POLYDRAW No DATA    
  {
    for (P in seq(1,length(PLCT3[,1]))){
      highlight.area=as.vector(PLCT3[P,c(1,2,2)])
      #print(highlight.area)
      start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[2]));
      end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[3]));
      
      polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 1]*innerPos);
      polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 2]*innerPos);
      
      
      PolyX<-rbind(PolyX,c(polygon.x))
      PolyY<-rbind(PolyY,c(polygon.y))
      #print(N)
      N<-N+1
    }
    
    PolyX<-PolyX[-1,]
    PolyY<-PolyY[-1,]
    PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
    PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
    print('yes')
    polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
    PaulN<-PaulN+0.05
  }  
  
}

#F1 - SNP
{
  PLCT<-read.table("/Volumes/UNTITLED/IMAD_2020/KDH/KDH_F1_P1_informative.vcf")
  PX<-c()
  ssC<-PLCT[,1]
  PLCT3<-c(0,0)
  for (i in unique(PLCT[,1])){
    PLCT1<-PLCT[PLCT[,1]==i,3]
    HP<-hist(PLCT1,plot = F,breaks = 250 )
    length(HP$counts)
    length(HP$breaks)
    PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
    PLCT3<-rbind(PLCT3,PLCT2)
  }
  PLCT3<-PLCT3[-1,]
  
  color.gradient <- function(x, colors=c("white","purple"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  CG<-color.gradient(HP$counts)
  #plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')
  
  
  setwd('/Volumes/UNTITLED/IMAD_2020/')
  #POLYDRAW SNP
  
  N<-1
  PaulN<-1.15
  print(N)
  PolyX<-c(0,0)
  PolyY<-c(0,0)
  innerPos<-PaulN
  outerPos<-PaulN + 0.03
  
  
  #POLYDRAW No DATA    
  {
    for (P in seq(1,length(PLCT3[,1]))){
      highlight.area=as.vector(PLCT3[P,c(1,2,2)])
      #print(highlight.area)
      start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                      as.numeric(highlight.area[2]));
      end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[3]));
      
      polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 1]*innerPos);
      polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                     RCircos.Pos[end.pos:start.pos, 2]*innerPos);
      
      
      PolyX<-rbind(PolyX,c(polygon.x))
      PolyY<-rbind(PolyY,c(polygon.y))
      #print(N)
      N<-N+1
    }
    
    PolyX<-PolyX[-1,]
    PolyY<-PolyY[-1,]
    PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
    PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
    print('yes')
    polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
    PaulN<-PaulN+0.05
  }  
  
}

#Beetle
{
Beetle7Den1<-read.table('/Volumes/UNTITLED/Projects/GSS/R_Circos/Tony_CHmarks/Paul_Blast/Beetle7.fa.out')
Beetle7Den1[,2]
ssC<-c()
for (ss in Beetle7Den1[,2]){
  ss<-strsplit((ss),'_')[[1]][1]
  ssC<-c(ssC,ss)
}
Beetle7Den1<-cbind(Beetle7Den1,ssC)
PLCT<-Beetle7Den1

PLCT3<-c(0,0)
for (i in unique(PLCT[,13])){
  PLCT1<-PLCT[PLCT[,13]==i,9]
  HP<-hist(PLCT1,plot = F,breaks = 25 )
  length(HP$counts)
  length(HP$breaks)
  PLCT2<-cbind(rep(i,length(HP$counts)),HP$breaks[seq(0,length(HP$breaks)-1)])
  PLCT3<-rbind(PLCT3,PLCT2)
}
PLCT3<-PLCT3[-1,]
PLCT3

color.gradient <- function(x, colors=c("white","red"), colsteps=1000) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
CG<-color.gradient(HP$counts)
#plot(LOS1$x,rep(0.9,250),pch=19,col='white',cex = 0.25,ylim=c(-1,1),main='ABBA BABA')


setwd('/Volumes/UNTITLED/IMAD_2020/')
#POLYDRAW SNP

N<-1
PaulN<-0.27
print(N)
PolyX<-c(0,0)
PolyY<-c(0,0)
innerPos<-PaulN
outerPos<-PaulN + 0.03


#POLYDRAW No DATA    
{
  for (P in seq(1,length(PLCT3[,1]))){
    highlight.area=as.vector(PLCT3[P,c(1,2,2)])
    #print(highlight.area)
    start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                    as.numeric(highlight.area[2]));
    end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                  as.numeric(highlight.area[2]));
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    
    PolyX<-rbind(PolyX,c(polygon.x))
    PolyY<-rbind(PolyY,c(polygon.y))
    #print(N)
    N<-N+1
  }
  
  PolyX<-PolyX[-1,]
  PolyY<-PolyY[-1,]
  PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
  PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
  print('yes')
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border=CG)
  PaulN<-PaulN+0.05
}  
}

#GeneDen
{
PolyX<-c(0,0)
PolyY<-c(0,0)
GeneDen<-read.table('/Volumes/UNTITLED/Projects/GSS/R_Circos/genes_total_circos_1.txt')
GeneDen
GeneDenTt<-c()
SSC<-c()
for (i in unique(GeneDen[,1])){
  GeneDen1<-hist(GeneDen[GeneDen[,1]==i,2],plot=F,breaks=100)
  GeneDen1$breaks<-GeneDen1$breaks[-length(GeneDen1$breaks)]
  GeneDenTt<-rbind(GeneDenTt,cbind(as.numeric(GeneDen1$breaks),as.numeric(GeneDen1$counts)))
  SSC<-c(SSC,rep(i,length(GeneDen1$breaks)))
}
noquote(SSC)
GeneDenTt<-cbind(SSC,GeneDenTt)
GeneDenTt
as.numeric(GeneDenTt[,2])

{  
  #HIST
  histValues <- as.numeric(GeneDenTt[,3]);
  max.value <- max(histValues);
  min.value <- min(histValues);
  histHeight <- histValues/max(histValues)
  histHeight
  
  histHeight
  GeneDenTt<-cbind(GeneDenTt,histHeight)
  GeneDenTt
  GeneDenTt[,c(1,2,2)]
  head(GeneDenTt)
  length(GeneDenTt[,1])
  
  
  
  for (P in seq(1,length(GeneDenTt[,1]))){
    OP<-0.4
    IP<-0.3
    IPOP<-as.numeric((OP-IP))
    innerPos<-IP
    outerPos<-as.numeric(IP)+(IPOP*as.numeric(GeneDenTt[P,4]))
    start.pos <- RCircos.Data.Point(as.character(GeneDenTt[P,1]),as.numeric(GeneDenTt[P,2]))
    end.pos <- RCircos.Data.Point(as.character(GeneDenTt[P,1]),as.numeric(GeneDenTt[P,2]));
    #print(start.pos)
    #print(end.pos)
    #print(innerPos)
    print(outerPos)
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    polygon(x=polygon.x, y=polygon.y,border = 'lightblue')
    
    polygon.x
    polygon.y
    PolyX<-rbind(PolyX,c(polygon.x))
    PolyY<-rbind(PolyY,c(polygon.y))
    
  }
  
  PolyX<-PolyX[-1,]
  PolyY<-PolyY[-1,]
  PX<-cbind(PolyX[,1],PolyX[,1],PolyX[,2],PolyX[,2])
  PY<-cbind(PolyY[,1],PolyY[,1],PolyY[,2],PolyY[,2])
  print('yes')
  polygon(x=cuts(t(PX)), y=cuts(t(PY)), border='lightblue')
  

}

}

#Repeats
{
RepetativeTt<-c(0,0,0)
Repetative<-read.table('/Volumes/UNTITLED/PacBio_EL_10/repeat_den_hist_1')
for (i in unique(Repetative[,1])){
  CentChr<-i
  Repetative1<-hist(Repetative[Repetative[,1]==i,2],plot=F,breaks=100)
  Repetative1$breaks<-Repetative1$breaks[-length(Repetative1)]
  RepetativeTt<-rbind(RepetativeTt,cbind(rep(i,length(Repetative1$breaks)),Repetative1$breaks,Repetative1$counts))  
}

RepetativeTt<-RepetativeTt[-1,]
RepetativeTt
#plot(Repetative1$breaks,Repetative1$counts,pch=19)

RepeatDenTT1<-data.frame()
RepeatDen<-read.table('/Volumes/UNTITLED/Projects/GSS/R_Circos/Tony_CHmarks/Paul_Blast/Beta_vulgaris_main_satellites.fa.out')
unique(RepeatDen[,1])
RepeatDen
ssC<-c()
for (ss in RepeatDen[,2]){
  ss<-strsplit((ss),'_')[[1]][1]
  ssC<-c(ssC,ss)
}
RepeatDen<-cbind(RepeatDen,ssC)
RepeatDen

for (N in unique(ssC)){
  for (i in unique(RepeatDen[,1])){
    RepeatDen1<-hist(RepeatDen[RepeatDen[,1]==i & RepeatDen[,13]==N,9],plot=F,breaks = 100)
    RepeatDen1$breaks<-RepeatDen1$breaks[-length(RepeatDen1$breaks)]
    for (xN in seq(1,length(RepeatDen1$breaks))){
      if (RepeatDen1$breaks[xN] < 50000){
        RepeatDenTT1<-rbind(RepeatDenTT1,cbind(i,as.numeric(RepeatDen1$breaks),
                                               as.numeric(RepeatDen1$breaks+50000),as.numeric(RepeatDen1$counts),N))
      }
      if (RepeatDen1$breaks[xN] >= 50000){
        RepeatDenTT1<-rbind(RepeatDenTT1,cbind(i,as.numeric(RepeatDen1$breaks)-50000,as.numeric(RepeatDen1$breaks),
                                               as.numeric(RepeatDen1$counts),N))
      }
    }
  }
}
}
RepeatDenTT1
RepeatDenTt



SSC<-c()
ShuJunTt<-data.frame()
ShuJun<-read.table('/Volumes/UNTITLED/Projects/GSS/EL10_fst/Global/ShuJun/Copia.gff')
ShuJun
for (i in unique(ShuJun[,1])){
  ShuJun1<-hist(ShuJun[ShuJun[,1]==i,4],plot=F,breaks=100)
  ShuJun1$breaks<-ShuJun1$breaks[-length(ShuJun1$breaks)]
  ShuJunTt<-rbind(ShuJunTt,cbind(as.numeric(ShuJun1$breaks),as.numeric(ShuJun1$counts)))
  SSC<-c(SSC,rep(i,length(ShuJun1$breaks)))
}
head(ShuJunTt)
ShuJunTt[,1]
noquote(SSC)
ShuJunTtCopia<-(cbind(SSC,ShuJunTt))


#Copia
{  
  #HIST
  histValues <- as.numeric(ShuJunTtCopia[,3]);
  max.value <- max(histValues);
  min.value <- min(histValues);
  histHeight <- ShuJunTtCopia[,3]/max(ShuJunTtCopia[,3])
  histHeight
  
  histHeight
  ShuJunTtCopia<-cbind(ShuJunTtCopia,histHeight)
  ShuJunTtCopia
  
  HistX<-data.frame()
  HistY<-data.frame()
  ShuJunTtCopia[,c(1,2,2)]
  head(ShuJunTtCopia)
  length(ShuJunTtCopia[,1])
  as.numeric(ShuJunTtCopia[,2])
  
  
  
  for (P in seq(1,length(ShuJunTtCopia[,1]))){
    OP<-0.55
    IP<-0.45
    IPOP<-as.numeric((OP-IP))
    innerPos<-IP
    outerPos<-as.numeric(IP)+(IPOP*as.numeric(ShuJunTtCopia[P,4]))
    start.pos <- RCircos.Data.Point(as.character(ShuJunTtCopia[P,1]),as.numeric(ShuJunTtCopia[P,2]))
    end.pos <- RCircos.Data.Point(as.character(ShuJunTtCopia[P,1]),as.numeric(ShuJunTtCopia[P,2]));
    #print(start.pos)
    #print(end.pos)
    #print(innerPos)
    print(outerPos)
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    
    polygon.x
    polygon.y
    
    HistX<-rbind(HistX,polygon.x)
    HistY<-rbind(HistX,polygon.y)
    polygon(polygon.x,polygon.y,border=adjustcolor('blue', alpha.f = 0.25))
    
  }
}





SSC<-c()
ShuJunTt<-data.frame()
ShuJun<-read.table('/Volumes/UNTITLED/Projects/GSS/EL10_fst/Global/ShuJun/Gypsy.gff')
ShuJun
for (i in unique(ShuJun[,1])){
  ShuJun1<-hist(ShuJun[ShuJun[,1]==i,4],plot=F,breaks=100)
  ShuJun1$breaks<-ShuJun1$breaks[-length(ShuJun1$breaks)]
  ShuJunTt<-rbind(ShuJunTt,cbind(as.numeric(ShuJun1$breaks),as.numeric(ShuJun1$counts)))
  SSC<-c(SSC,rep(i,length(ShuJun1$breaks)))
}
head(ShuJunTt)
ShuJunTt[,1]
noquote(SSC)
ShuJunTtGypsy<-as.data.frame(cbind(SSC,ShuJunTt))

#Gypsy
{
  #HIST
  histValues <- as.numeric(ShuJunTtGypsy[,3]);
  max.value <- max(histValues);
  min.value <- min(histValues);
  histHeight <- ShuJunTtGypsy[,3]/ShuJunTtGypsy[order(ShuJunTtGypsy[,3],decreasing=T)[3],3]
  histHeight
  
  histHeight
  ShuJunTtGypsy<-cbind(ShuJunTtGypsy,histHeight)
  ShuJunTtGypsy
  
  HistX<-data.frame()
  HistY<-data.frame()
  ShuJunTtGypsy[,c(1,2,2)]
  head(ShuJunTtGypsy)
  length(ShuJunTtGypsy[,1])
  as.numeric(ShuJunTtGypsy$V2)
  
  
  
  for (P in rownames(ShuJunTtGypsy)){
    OP<-0.55
    IP<-0.45
    IPOP<-as.numeric((OP-IP))
    innerPos<-IP
    outerPos<-as.numeric(IP)+(IPOP*as.numeric(ShuJunTtGypsy[P,4]))
    start.pos <- RCircos.Data.Point(as.character(ShuJunTtGypsy[P,1]),as.numeric(ShuJunTtGypsy[P,2]))
    end.pos <- RCircos.Data.Point(as.character(ShuJunTtGypsy[P,1]),as.numeric(ShuJunTtGypsy[P,2]));
    #print(start.pos)
    #print(end.pos)
    #print(innerPos)
    print(outerPos)
    
    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                   RCircos.Pos[end.pos:start.pos, 2]*innerPos);
    
    
    polygon.x
    polygon.y
    
    HistX<-rbind(HistX,polygon.x)
    HistY<-rbind(HistX,polygon.y)
    polygon(polygon.x,polygon.y,col='blue',border=adjustcolor('red', alpha.f = 0.25))
    
  }
}








#plot(Beetle7Den$breaks,Beetle7Den$counts,pch=19)
##PLOT##

dev.copy(png,"/Volumes/UNTITLED/IMAD_2020/KDH/KDH13.png",height=11,width=11,units='in',res=1000)
dev.off()
