library('survival')
library('survminer')
library(dplyr)

# Load Gene Expression cohort data
# Load the Gene expression data. Each row corresponds to gene while each column corresponds to patients.
data= read.table("GeneExpression_log2Counts_TNBC_HER2pos_ERpos_Normal_normalized.txt", sep="\t", header=T,check.names=F, stringsAsFactor=F)

# Load sample information data
# This file contains information of each patient, 
sampleinfo = read.table("sample_group_info.txt", sep="\t", header=T,check.names=F, stringsAsFactor=F)
rownames(sampleinfo) = sampleinfo$sampleName;
sampleinfo$patientid = substr(sampleinfo$sampleName,1, 12);

# Stage summary
table(sampleinfo$stageCode)

# EPHsubtype
table(sampleinfo$groupName)

# Load Clinical Information of each patient
clinical = read.table("TCGA-BRCA_clinical.csv", sep=",", header=T,check.names=F, stringsAsFactor=F)
rownames(clinical) = clinical$submitter_id;

# Map TCGA patiend id to EPH status (ERpos/TNBC/HER2pos/Normal)
ephstatus = sampleinfo[colnames(data)[2:ncol(data)],"groupName"];
names(ephstatus) = colnames(data)[2:ncol(data)];


# First gene or 1st row
rowno=1;

# All sample names
samplename = colnames(data)[2:ncol(data)]

# All patient names
patientname = sampleinfo[samplename,"patientid"];

# Map clinical data with gene expression data
tempclindata = clinical[,c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")];
tempclindata = tempclindata[patientname,];
tempclindata$ephstatus = ephstatus[samplename];
tempclindata$genexp = unlist(data[rowno,samplename]);

notDead <- is.na(tempclindata$days_to_death)

if (any(notDead == TRUE)) {
       tempclindata$days_to_death[notDead] <- tempclindata[notDead, "days_to_last_follow_up"]
}

tempclindata$s <- grepl("dead", tempclindata$vital_status, ignore.case = TRUE)
tempclindata$ephstatus <- as.factor(tempclindata$ephstatus)

## Fit Proportional Hazards Regression Model
res.cox <- coxph(Surv(days_to_death,s) ~ ephstatus + genexp, data =  tempclindata)
summary(res.cox)

## Do processing
newdata=data.frame(ephstatus=unique(as.character(tempclindata$ephstatus)), genexp = sapply(unique(as.character(tempclindata$ephstatus)), function(x) { filter(tempclindata, ephstatus==x) %>% select(genexp) %>% unlist %>% mean;}));

## Create survival curves
fit <- survfit(res.cox, newdata = newdata)

## Drawing Survival Curves Using ggplot2
ggsurvplot(fit, data = newdata, conf.int = TRUE, legend.labs=newdata$ephstatus,ggtheme = theme_minimal())
        
           


   
   
   
        
        
    
    




