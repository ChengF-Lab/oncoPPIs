uniprot_length <- read.table("UniProt_len.txt",colClasses=c("character","numeric"))

psites = list()
data <- read.table('PPI_sites.txt',sep="\t",colClasses="character")  
for(i in 1:nrow(data)){
    psites[[data[i,1]]] <- unlist(strsplit(data[i,2],","))
}

cancer_types <- sapply(list.files('HQinterfaces_mutations/'),function(x) unlist(strsplit(x,split="_"))[1])

for(cancer in cancer_types){
    Mutations <- read.table(paste0('/data/Data/TCGA_Somatic_Mutations/Mutect2/TCGA.', cancer, '.mutect.NonSilent.maf'),
        skip=1,header=TRUE,colClasses="character",sep="\t")
    PPI_mutations <- read.table(paste0('HQinterfaces_mutations/', cancer, '_HQinterfaces_mutations.txt'),colClasses="character")
    all_uniprot <- unique(PPI_mutations[,3])
    uniprot_gene_mapping <- rep(NA,length(all_uniprot))
    names(uniprot_gene_mapping) <- all_uniprot
    uniprot_pvalues <- rep(1,length(all_uniprot))
    names(uniprot_pvalues) <- all_uniprot
    for(uniprot in all_uniprot){
        uniprot_gene_mapping[uniprot] <- PPI_mutations[match(uniprot,PPI_mutations[,3]),2]
        mutation_number <- sum(Mutations[,68]==uniprot)
        PPI_mutation_number <- sum(PPI_mutations[,3]==uniprot)
        gene_len <- uniprot_length[match(uniprot,uniprot_length[,1]),2]
        bind_len <- length(psites[[uniprot]])
        if(PPI_mutation_number > 0 & mutation_number > 0){
            p <- binom.test(x=PPI_mutation_number,n=mutation_number,p=bind_len/gene_len,alternative="greater")$p.value
            uniprot_pvalues[uniprot] <- p
        }
    }
    write.table(cbind(uniprot_gene_mapping,uniprot_pvalues),
                file=paste("./uniprot_pvalues/",cancer,"_uniprot_pvalues.txt",sep=""),
                col.names=FALSE,
                sep="\t",
                quote=FALSE)
}

cancer <- cancer_types[1]
Mutations <- read.table(paste0('/data/Data/TCGA_Somatic_Mutations/Mutect2/TCGA.', cancer, '.mutect.NonSilent.maf'),
        skip=1,header=TRUE,colClasses="character",sep="\t")
PPI_mutations <- read.table(paste0('HQinterfaces_mutations/', cancer, '_HQinterfaces_mutations.txt'),colClasses="character")
#Pan-cancer
for(cancer in cancer_types[2:33]){
    Mutations <- rbind(Mutations,
                       read.table(paste0('/data/Data/TCGA_Somatic_Mutations/Mutect2/TCGA.', cancer, '.mutect.NonSilent.maf'),
                       skip=1,header=TRUE,colClasses="character",sep="\t"))
    PPI_mutations <- rbind(PPI_mutations,
                           read.table(paste0('HQinterfaces_mutations/', cancer, '_HQinterfaces_mutations.txt'),
                            colClasses="character"))
}
all_uniprot <- unique(PPI_mutations[,3])
uniprot_gene_mapping <- rep(NA,length(all_uniprot))
names(uniprot_gene_mapping) <- all_uniprot
uniprot_pvalues <- rep(1,length(all_uniprot))
names(uniprot_pvalues) <- all_uniprot
for(uniprot in all_uniprot){
    uniprot_gene_mapping[uniprot] <- PPI_mutations[match(uniprot,PPI_mutations[,3]),2]
    mutation_number <- sum(Mutations[,68]==uniprot)
    PPI_mutation_number <- sum(PPI_mutations[,3]==uniprot)
    gene_len <- uniprot_length[match(uniprot,uniprot_length[,1]),2]
    bind_len <- length(psites[[uniprot]])
    if(PPI_mutation_number > 0 & mutation_number > 0){
        p <- binom.test(x=PPI_mutation_number,n=mutation_number,p=bind_len/gene_len,alternative="greater")$p.value
        uniprot_pvalues[uniprot] <- p
    }
}
write.table(cbind(uniprot_gene_mapping,uniprot_pvalues),
            file=paste("./uniprot_pvalues/Pan-cancer_uniprot_pvalues.txt",sep=""),
            col.names=FALSE,
            sep="\t",
            quote=FALSE)
