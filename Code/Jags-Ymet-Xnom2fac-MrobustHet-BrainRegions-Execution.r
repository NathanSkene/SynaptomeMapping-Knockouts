# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
a.c <- function(x){return(as.character(x))}


##########################################################################################
# SET PARAMETERS:
##########################################################################################
# << THE BELOW LINE WILL NEED EDITING FOR ANYONE ELSE RUNNING THIS CODE ON WINDOWS >>
baseDir = "~/ownCloud/Shared with Nathan/Synaptome Map Statistics/Knockout/"
setwd(baseDir)
# << EDIT BELOW LINE TO CONTAIN NAMES OF THE DATA FILES >>
file_tag = c("Density_GHB","Density_GGW")
graphFileType = "pdf" 

#DEBUGGING=TRUE
DEBUGGING=FALSE
load_ref = "Density_GHB-Density_GGW_2016-04-07 13_45_06"

##########################################################################################
# SETUP FOLDERS TO STORE THE RESULTING DATA
##########################################################################################
folder_ref = gsub(":","_",sprintf("%s_%s",paste(file_tag,collapse="-"),Sys.time()))
analysis_dir = sprintf("%sOutput/%s",baseDir,folder_ref)
if(DEBUGGING==TRUE){
	analysis_dir = sprintf("%sOutput/%s",baseDir,load_ref)	
}
dir.create(file.path(sprintf("%sOutput",baseDir), folder_ref))
dir.create(file.path(analysis_dir, "Figs_DiagMCMC"))
dir.create(file.path(analysis_dir, "Figs_PlotMCMC"))
dir.create(file.path(analysis_dir, "Figs_RegionPlots"))
dir.create(file.path(analysis_dir, "mcmcChains"))
dir.create(file.path(analysis_dir, "ResultsTables"))

##########################################################################################
# STEP 1: LOAD AND PREPARE THE DATA
##########################################################################################

FULL_RES = data.frame()
library(reshape)
variable_type = gsub("_.*","",file_tag)[1]
synapse_data=data.frame()
colonies = gsub(".*_","",file_tag)
for(ff in file_tag){
	colony = gsub(".*_","",ff)
	raw_synapse_data = read.csv(sprintf("Data/%s.csv",ff))
	raw_synapse_data = t(raw_synapse_data)
	colnames(raw_synapse_data)=raw_synapse_data[1,]
	raw_synapse_data = raw_synapse_data[-1,]
	raw_synapse_data = cbind(mID=rownames(raw_synapse_data),raw_synapse_data)
	raw_synapse_data = raw_synapse_data[,-which(colnames(raw_synapse_data)=="Parents pair")]
	syn_dat = melt(data.frame(raw_synapse_data),id.vars=c("mID","SEX","GENOTYPE","COLONY","Litter"))
	syn_dat$value = as.numeric(as.character(syn_dat$value))
	syn_dat$region = gsub("_.*","",syn_dat$variable)
	syn_dat = syn_dat[!is.na(syn_dat$value),]
	syn_dat$value = as.numeric(syn_dat$value)
	syn_dat$GENOTYPE = as.factor(gsub("het/","",a.c(syn_dat$GENOTYPE)))
	syn_dat$GENOTYPE = sprintf("%s %s",colony,syn_dat$GENOTYPE)
	if(dim(synapse_data)[2]==0){
		synapse_data = syn_dat
	}else{
		synapse_data = rbind(synapse_data,syn_dat)
	}
}

########################################################################################################
# STEP 2: FOR EACH FULL REGION (e.g. hippocampus/cortex/striatum), SEQUENTIALLY PERFORM THE ANALYSIS
# ---------- NOTE, SUBREGIONS, e.g. Hippocampus_DGgrSup ARE ANALYSED WITHIN THEIR ASSOCIATED FULL REGION
########################################################################################################

                                            
for(region in unique(synapse_data$region)){ 
	######################################################################################
	### SEPERATE THE DATA RELEVANT TO THAT REGION ########################################
	######################################################################################
    myDataFrame = synapse_data[synapse_data$region==region,]
    myDataFrame$variable = as.factor(as.character(myDataFrame$variable))
    myDataFrame$GENOTYPE = as.factor(myDataFrame$GENOTYPE)
    
    # Specify the column names in the data file relevant to the analysis:
    x1Name="GENOTYPE" 
    x2Name="variable" 
    yName="value" 
    
    # Specify contrasts
    contrasts = list() 
    x1contrasts = list()
    for(cc in colonies){
    	x1contrasts[[length(x1contrasts)+1]] = list( sprintf("%s wt",cc) , sprintf("%s hom",cc) , compVal=0.0 , ROPE=c(-0.1,0.1) )
    }
    
    
    ##########################################################################################
    # LOAD AND RUN THE BAYESIAN MODELLING FOR THIS REGION
    ##########################################################################################

    #------------------------------------------------------------------------------- 
    # Load the relevant model into R's working memory:
    source("Code/Jags-Ymet-Xnom2fac-MrobustHet-BrainRegions-Model.r")
    #------------------------------------------------------------------------------- 
    # Generate the MCMC chain:
    if(DEBUGGING==FALSE){
	    mcmc_fileNameRoot = sprintf("%s/mcmcChains/mcmc_%s_%s_",analysis_dir,region,variable_type) 
	    mcmcCoda = genMCMC( datFrm=myDataFrame , 
	                        yName=yName , x1Name=x1Name , x2Name=x2Name ,
	                        numSavedSteps=15000 , thinSteps=5 , saveName=mcmc_fileNameRoot )
	                        # SET numSavedSteps=15000
	    save(mcmcCoda,file=sprintf("%s/mcmcChains/mcmcCoda_%s_%s.RData",analysis_dir,variable_type,region))
    }else{
    	load(file=sprintf("%s/mcmcChains/mcmcCoda_%s_%s.RData",analysis_dir,variable_type,region))	
    }     
    
    ##########################################################################################
    # COMMENTED SECTIONS ARE FOR CHECKING THE MODELLING OF PARAMETERS WORKED WELL
    ##########################################################################################
    #------------------------------------------------------------------------------- 
    DiagMCMC_fileNameRoot = sprintf("%s/Figs_DiagMCMC/%s_%s",analysis_dir,region,variable_type)
    ## Display diagnostics of chain, for specified parameters:
    parameterNames = varnames(mcmcCoda) 
    show( parameterNames ) # show all parameter names, for reference
    #for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma[1,1]","ySigma[1,7]","ySigma[5,7]","nu") ) {
    for ( parName in parameterNames) {
        diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=DiagMCMC_fileNameRoot , saveType=graphFileType )
    }
    #------------------------------------------------------------------------------- 
    PlotMCMC_fileNameRoot = sprintf("%s/Figs_PlotMCMC/%s_%s",analysis_dir,region,variable_type)
    ## Get summary statistics of chain:
	#     summaryInfo = smryMCMC( mcmcCoda , 
	#                             datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name ,
	#                             saveName=fileNameRoot )
	#     show(summaryInfo)
	#     # Display posterior information:
	#     plotMCMC( mcmcCoda , 
	#               datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name ,
	#               saveName=PlotMCMC_fileNameRoot , saveType=graphFileType )
    #------------------------------------------------------------------------------- 
    
    
    ##########################################################################################
    # FOR EACH KNOCKOUT LINE...
    ##########################################################################################
    # Other specific comparisons of cells:
    region_plot_base = sprintf("%s/Figs_RegionPlots/RegionPlot_%s",analysis_dir,variable_type)
    for(cc in colonies){
   	    ##########################################################################################
	    # HDI & ROPE ANALYSIS FOR EACH SUBREGION WITHIN THIS REGION+KNOCKOUT
	    ##########################################################################################

        print("FIRST, SUBREGIONS")
    	regions = unique(myDataFrame$variable)
        THISx1 = sprintf("%s hom",cc)
        THATx1 = sprintf("%s wt",cc)	
        THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
        THATidx = which(levels(myDataFrame[,x1Name])==THATx1)    
    	for(reg in regions){
   		  # Get Scale (Standard Deviation) of the T-Distribution for that brain region
   		  # - Remember, scale is only calculated on the level of region, so is equal for the two genotypes	
    	  ATx2 = reg
    	  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
    	  if(length(regions)>1){
    	  	sd_dist = as.matrix(mcmcCoda)[,paste("ySigma[",ATidx,"]",sep="")]
    	    reg_sd  = median(as.matrix(mcmcCoda)[,paste("ySigma[",ATidx,"]",sep="")])/3
    	  }else{
    	      reg_sd  = median(as.matrix(mcmcCoda)[,"ySigma"])/3  
    	  }
		  estimate_mode <- function(x) {
			  d <- density(x)
			  d$x[which.max(d$y)]
		  }      	  
		  scale_stats = data.frame(scale_mean=mean(sd_dist),scale_median=median(sd_dist),scale_mode=estimate_mode(sd_dist))

    	  # Get (and plot) the difference between the two means
    	  # THIS x1level minus THAT x1level at AT x2level:
    	  meanDiffs = as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
    		  as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] 
    	  pdf(sprintf("%s_%s_%s.pdf",region_plot_base,cc,reg),height=4,width=5)
    	  compInfo = plotPost( 
    		meanDiffs, 
    		main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    		xlab=paste("Difference in",yName) , 
    		compVal=0 ,ROPE=c(-reg_sd,reg_sd) )
    	  show(compInfo)
    	  dev.off()
    	  
    	  # Calculate the Cohen's D distribution
  	  
    	  cohensD = meanDiffs/sd_dist
    	  cohensD_quantiles = quantile(cohensD,c(0.05,0.5,0.95))
    	  cohensD_stats = data.frame(cohensD_mode=estimate_mode(cohensD),cohensD_mean=mean(cohensD),cohensD_point05quantile=cohensD_quantiles[1],cohensD_median=cohensD_quantiles[2],cohensD_point95quantile=cohensD_quantiles[3])
    	  
    	  add_dat = cbind(region=ATx2,regionType="SubRegion",colony=cc,compInfo,cohensD_stats,scale_stats)
    	  if(dim(FULL_RES)[1]==0){
    		FULL_RES = add_dat
    	  }else{
    		FULL_RES = rbind(FULL_RES,add_dat)
    	  }
    	}
    	
   	    ##########################################################################################
	    # HDI & ROPE ANALYSIS FOR EACH GROUP OF SUBREGIONS (e.g. CA3) WITHIN THIS REGION+KNOCKOUT
	    ##########################################################################################    	
    	#------------------------------------------------------------------------------- 
    	# Look at overall regions
        print("THEN, GROUPED REGIONS")
    	subregs=c("Overall")
    	if(region=="Hippocampus"){subregs=c(subregs,"CA1","CA2","CA3","DG")}
        if(region=="CorticalSubplate"){subregs=c(subregs,"BLA","EP")}
        if(region=="Neocortex"){subregs=c(subregs,"AUD","MO","RSP","SS")}
        if(region=="OlfactoryAreas"){subregs=c(subregs,"COA","PAA","PIR")}
                             
    	for(subR in subregs){
   		  # Get Scale (Standard Deviation) of the T-Distribution for that brain region
   		  # - Remember, scale is only calculated on the level of region, so is equal for the two genotypes	
    	  ATx2 = sprintf("%s %s",region,subR)
    	  if(subR!="Overall"){ATidx   = grep(subR,levels(myDataFrame[,x2Name]))}
    	  if(subR=="Overall"){ATidx   = 1:length(levels(myDataFrame[,x2Name]))}
    	  lenSR   = length(ATidx)  
    	  if(length(regions)>1){
    	  	sd_dist = as.matrix(mcmcCoda)[,paste("ySigma[",ATidx,"]",sep="")]
    	    reg_sd  = median(apply(as.matrix(mcmcCoda)[,paste("ySigma[",ATidx,"]",sep="")],1,mean))/3
    	  }else{
    	  	sd_dist = as.matrix(mcmcCoda)[,"ySigma"]
    	    reg_sd  = median(as.matrix(mcmcCoda)[,"ySigma"])/3
    	  }
		  scale_stats = data.frame(scale_mean=mean(sd_dist),scale_median=median(sd_dist),scale_mode=estimate_mode(sd_dist))
    	  
    	  # Get (and plot) the difference between the two means
    	  # THIS x1level minus THAT x1level at AT x2level:    	  
    	  pdf(sprintf("%s_%s_%s.pdf",region_plot_base,cc,reg),height=4,width=5)
    	  avg_reg_wt  = as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx[1],"]",sep="")]
    	  avg_reg_hom = as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx[1],"]",sep="")]	  
    	  if(lenSR>1){
        	  for(i in 2:lenSR){	
        		  avg_reg_wt  = avg_reg_wt + as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx[i],"]",sep="")]
        		  avg_reg_hom = avg_reg_hom + as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx[i],"]",sep="")]	  
        	  }
    	  }
    	  avg_reg_wt = avg_reg_wt/lenSR
    	  avg_reg_hom = avg_reg_hom/lenSR	  
    	  compInfo = plotPost( 
    		avg_reg_wt - avg_reg_hom , 
    		main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    		xlab=paste("Difference in",yName) , 
    		compVal=0 ,ROPE=c(-reg_sd,reg_sd) )
    	  show(compInfo)
    	  dev.off()
    	  meanDiffs = avg_reg_wt - avg_reg_hom
    	  
    	  # Calculate the Cohen's D distribution
    	  cohensD = meanDiffs/sd_dist
    	  cohensD_quantiles = quantile(cohensD,c(0.05,0.5,0.95))
    	  cohensD_stats = data.frame(cohensD_mode=estimate_mode(cohensD),cohensD_mean=mean(cohensD),cohensD_point05quantile=cohensD_quantiles[1],cohensD_median=cohensD_quantiles[2],cohensD_point95quantile=cohensD_quantiles[3])
    	  
    	  
    	  # Add rows to FULL_RES
    	  add_dat = cbind(region=ATx2,regionType="GreaterRegion",colony=cc,compInfo,cohensD_stats,scale_stats)
    	  FULL_RES = rbind(FULL_RES,add_dat)
    	}
    }
}
##########################################################################################
# PREPARE THE DATA FOR SAVING
##########################################################################################
FULL_RES = data.frame(FULL_RES)
a.n <- function(x){as.numeric(a.c(x))}
q = p.adjust(c(1-a.n(FULL_RES$pLtROPE),1-a.n(FULL_RES$pGtROPE)),method="BH")
FULL_RES$qROPE_Ls = q[1:dim(FULL_RES)[1]]
FULL_RES$qROPE_Gt = q[(dim(FULL_RES)[1]+1):length(q)]
FULL_RES$minP = apply(cbind(1-a.n(FULL_RES$pLtROPE),1-a.n(FULL_RES$pGtROPE)),1,min)
FULL_RES$minQ = apply(FULL_RES[,c("qROPE_Ls","qROPE_Gt")],1,min)
write.csv(FULL_RES,file=sprintf("%s/ResultsTables/FULL_RES_%s.csv",analysis_dir,variable_type))