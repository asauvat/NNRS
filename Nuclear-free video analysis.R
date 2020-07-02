library(MetaxpR)
library(MorphR)
library(EBImage)
library(RBioFormats)
#
library(doParallel)
library(pbapply)
#
library(magrittr)
library(RNiftyReg)

#=================================================================================================================================================================#
SERVER = 'VIDEO-SERVER'
DB = GetMDCInfo(SERVER)
#
PID = 587
PLT = DB[which(DB$PlateID==PID),]
PIF = GetPINFO(SERVER,PlateID=PID)
alls = lapply(1:nrow(PIF$SID),function(i)PIF$SID[i,])
#
INF = lapply(GetIMPath(SERVER,PlateID=PID,WellID=PIF$WID[1],SiteID=alls[[1]],TimePoint=1),function(x)read.metadata(x))
#
#IXMC-----
cn = sapply(strsplit(enc2utf8(sapply(INF,function(x)x$globalMetadata$`image-name`)),'[<]'),function(x)x[[1]]);names(cn)=names(INF)=cn
PS = eval(parse(text=INF[[1]]$globalMetadata$`spatial-calibration-x`))
#IXM XLS-
#cn = sapply(strsplit(enc2utf8(sapply(INF,function(x)x$seriesMetadata$Name)),'[<]'),function(x)x[[1]]);names(cn)=names(INF)=cn
#PS = eval(parse(text=INF[[1]]$seriesMetadata$XCalibration))

#=================================================================================================================================================================#

PM='/media/Hcs-screen10-vi/D/PlateMap'
DrugID = read.csv(file.path(PM,SERVER,paste0(PID,'.csv')),stringsAsFactors=F)

#=================================================================================================================================================================#
cp = 'TOR'
cl=makeCluster(15)
invisible({clusterEvalQ(cl,{library(MetaxpR);library(MorphR);library(EBImage)})})
clusterExport(cl,c('DrugID','cn','SERVER','PID','alls','cp'))
tstart=1;tstop=PIF$TID[2] #Skipping t0 because x/y shift
rgi = pblapply(tstart:tstop,function(ti){
  woi = subset(DrugID,grepl(cp,Drug))[,'WellID']
  IM = lapply(seq_along(cn),function(i)readImage(sapply(woi,function(w)GetIMPath(SERVER,PlateID=PID,WellID=w,SiteID = alls[[2]],TimePoint=ti)[i])))
  rg = lapply(IM,function(x){qu=quantile(x,probs = seq(0,1,5*10**-4));return(qu[c(1,length(qu)-1)])});names(rg)=cn
  return(rg)
},cl=cl);names(rgi)=tstart:tstop
stopCluster(cl);rm(cl)
#
save.image('PDATA.RData')

#=================================================================================================================================================================#

dir.create('SEGS')
sapply(tstart:tstop,function(ti)dir.create(file.path('SEGS',paste0('TimePoint_',ti)))) #Create folders

#=================================================================================================================================================================#
nCores=15

dat=list()  
for(ti in tstart:tstop){
  cl = makeCluster(nCores)
  clusterExport(cl,setdiff(ls(),'DB'))
  invisible({clusterEvalQ(cl,{library(MetaxpR);library(EBImage);library(MorphR);library(magrittr);library(S4Vectors)})})
  print(paste0('Analyzing TimePoint ',ti,'...'))
  tdat = pblapply(sort(PIF$WID),function(w){
    wdat = lapply(alls,function(sxy){
      
      #Read images ----------------------------------------------------------
      IP = GetIMPath(SERVER,PlateID=PID,WellID=w,SiteID=sxy,TimePoint=ti)
      IM = lapply(IP,function(x)suppressWarnings(readImage(x)));names(IM)=cn
      IN = lapply(cn,function(ci)normalize(IM[[ci]],inputRange=rgi[[as.character(ti)]][[ci]]))
        
      #CYTOPLASM-------------------------------------------------------------
      CG=IN$GFP-whiteTopHat(IN$GFP,makeBrush(25,'disc'));CG[which(CG<0)]=0
      #CG = CG-gblur(CG,300);CG[which(CG<0)]=0
      #
      SG = sigmoNormalize(CG,scaled=T,z0=4.5*10**-2,lambda=10**3)
      SM = SG>(otsu(SG)*0.95)
      SM = propagate(SM,opening(SM,makeBrush(51,'disc')),SM)
      #
      BTP = LowPass(blackTopHat(CG,makeBrush(5,'disc')))
      BTP = sigmoNormalize(BTP,scaled=F,z0=5.5*10**-1,lambda=10**3)
      BM = bwlabel(BTP>(otsu(BTP)*0.95))
      BM = rmObjects(BM,which(sapply(EBImage:::splitObjects(BM),length)<100))
      #
      CM = SM&!BM
      CM = watershed(distmap(CM),tolerance=5)
      #
      BD = dilate(Lap(CM)>1,makeBrush(3,'disc'))
      CM = bwlabel(fillHull(CM&!BD))
      #display(paintObjects(CM,toRGB(autoNormalize(IM$GFP))))
      
      #Dots segmentation-----------------------------------------------------
      ds = c(3,5,7,9,11,13,15)
      TP = Reduce(function(x,y)pmax(x,y),lapply(ds,function(x)whiteTopHat(IN$GFP,makeBrush(x,'disc'))))
      TP = sigmoNormalize(TP,scaled=T,z0=7.5*10**-2,lambda=10**4)
      #
      DM = TP>(otsu(TP)*0.9)
      DM = (DM&(erode(CM,makeBrush(15,'disc'))))*CM
      #
      DMp = bwlabel(DM>0)         
      
      #Final Image----------------------------------------------------------
      #rgbImage(green=Normalize(IM$GFP,autoRange=T,verbose=F))%>%paintObjects(x=CM,col=c('white',NA),thick=T)%>%paintObjects(x=DM,col=c('blue',NA))%>%display()
      
      #Export Seg----------------------------------------------------------
      SaveSeg(CM, ExportDir = 'SEGS',FileName=paste0('TimePoint_',ti,'/Cyto_',w,'_sx',sxy[1],'_sy',sxy[2]))
      SaveSeg(DM, ExportDir = 'SEGS',FileName=paste0('TimePoint_',ti,'/DotsG_',w,'_sx',sxy[1],'_sy',sxy[2]))
      
      #Features-------------------------------------------------------------
      
      FF = data.frame(PlateID=PID,Time=ti,WellID=w,SiteX=sxy[1],SiteY=sxy[2])
      
      CO = EBImage:::splitObjects(CM)
      if(length(CO)!=0){
        
        CF = data.frame(CellID=names(CO),do.call('rbind',lapply(CO,function(x){c(CytoArea=length(x)*(PS**2),GFPCytoIntensity=mean(IM$GFP[arrayInd(ind=x,.dim=dim(IM$GFP))]),
                                                                                 apply(arrayInd(ind=x,.dim=dim(IM$GFP)),2,mean))})))
        colnames(CF)[c(ncol(CF)-1,ncol(CF))]=c('PixelX','PixelY')
        
        #Dots: surface and count-------#
        DO = EBImage:::splitObjects(DM)
        DF = data.frame(CellID=names(DO),do.call('rbind',lapply(DO,function(x){c(DotsSurf=length(x)*(PS**2),DotsCount=length(unique(c(DMp[arrayInd(x,.dim=dim(IM$GFP))]))))})))
        
        #Dots: distance to their centroid#
        DOp = EBImage:::splitObjects(DMp)
        DDF = sapply(DO,function(x){
          cen = apply(arrayInd(ind=x,.dim=dim(IM$GFP)),2,mean) #Get dots centroid in cell
          sb = DOp[as.character(sort(unique(c(DMp[arrayInd(x,.dim=dim(IM$GFP))]))))] #Get labeled dots
          #
          d = sapply(sb,function(y){dist(rbind(cen,apply(arrayInd(ind=y,.dim=dim(IM$GFP)),2,mean)))})*PS #distance of dots from centroid
          ity = sapply(sb,function(y){sum(IM$GFP[y])}) #Integrated intensity of dots
          #
          pmean = sum(d*ity)/sum(ity) #Pondered mean of dots distance
          return(pmean)
        })
        DDF = data.frame(CellID=names(DDF),DD = DDF)
        
        #COMBINE-------------------------#
        FF = data.frame(PlateID=PID,Time=ti,WellID=w,SiteX=sxy[1],SiteY=sxy[2],Reduce(function(x,y)merge(x,y,by='CellID',all=T),list(CF,DF,DDF)))
        FF[,grep('Dot',colnames(FF))] = apply(FF[,grep('Dot',colnames(FF))],2,function(x){x[is.na(x)]=0;return(x)})
      }
      return(FF)
    })
    wdat=do.call(gtools::smartbind,wdat)
    return(wdat)
  },cl=cl)
  stopCluster(cl);rm(cl)
  tdat=do.call('rbind',tdat)
  dat=c(dat,list(tdat));rm(tdat)
};rm(ti)
save.image("SEGDATA.RData")
