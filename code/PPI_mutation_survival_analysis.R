library(RColorBrewer)
library(survival)

#cancer_types <- c("UCEC","SKCM","LUAD","COAD","BRCA")
#cancer_types <- sapply(list.files("HQinterfaces_mutation_edges/"),function(x) unlist(strsplit(x,split="_"))[1])
cancer_types <- c("COAD","LUSC")
for (cancer in cancer_types){
    clinical_file <- paste("~/Google Drive/Clinical/nationwidechildrens.org_clinical_patient_",tolower(cancer),".txt",sep="")
    if(file.exists(clinical_file)){
        clinical <- read.table(clinical_file,
                           skip=1,
                           sep="\t",
                           header=TRUE,
                           row.names=1,
                           quote="",
                           fill=TRUE)
        clinical <- clinical[-1,]
    
        pdf(paste0("survival_analysis/",cancer,"_PPI_mutation_survival_analysis.pdf"),width=5,height=5)
        sink(paste0("survival_analysis/",cancer,"_PPI_mutation_survival_analysis.txt"))
        mutated_PPI <- list()
        mutations <- read.table(paste0("HQinterfaces_mutation_edges/",cancer,"_HQinterfaces_mutation_edges.txt"),
        header=FALSE,sep="\t",fill=TRUE)

        for(i in 1:nrow(mutations)){
            Sample <- mutations[i,1]
            edges <- unlist(sapply(mutations[i,5],function(x) strsplit(x,',')))
            for(edge in edges){
                if(edge %in% names(mutated_PPI)){
                    mutated_PPI[[edge]] <- c(mutated_PPI[[edge]],Sample)
                }
                else{
                    mutated_PPI[[edge]] <- c(Sample)
                }
            }
        }

    days_to_followup <- apply(clinical,1,function(x) if(x["days_to_last_followup"]!="[Not Available]") as.numeric(x["days_to_last_followup"]) else as.numeric(x["days_to_death"]))
    vital_status <- apply(clinical,1,function(x) if(x["vital_status"]=="Dead") 1 else 0)
    
    for(feature in names(mutated_PPI)){
        if(length(mutated_PPI[[feature]])>=10){
            mutated_samples <- mutated_PPI[[feature]]
            mutated_samples <- intersect(rownames(clinical),mutated_samples)
            nonmutated_samples <- setdiff(rownames(clinical),mutated_samples)
            status <- c(rep("Mutated",length(mutated_samples)),rep("Wildtype",length(nonmutated_samples)))
            sdf <- survdiff( Surv(days_to_followup[c(mutated_samples,nonmutated_samples)], vital_status[c(mutated_samples,nonmutated_samples)])~ status)
            p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
            
            psite.surv <- survfit( Surv(days_to_followup[c(mutated_samples,nonmutated_samples)], vital_status[c(mutated_samples,nonmutated_samples)])~ status, conf.type="none")
            plot(psite.surv,  
                 xlab="Time", 
                 ylab="Survival Probability",
                 col=c("#ff7300","#64bbe3"),
                 bty='l',
                 lwd=3,
                 cex=1.5,
                 cex.axis=1,
                 cex.lab=1.5)
            title(main=paste(feature,paste("p-value=",round(p.val,2)),sep=":"))
            cat(feature)
            cat("\t")
            cat(p.val)
            cat("\n")
        }
    }
    dev.off()
    sink()
}
}

