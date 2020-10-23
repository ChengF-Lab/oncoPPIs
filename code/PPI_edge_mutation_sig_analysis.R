protein_coding_genes <- as.matrix(read.table("../snp_pruning/protein-coding_gene.txt",
                                   header=TRUE,sep="\t",quote="",fill=TRUE,colClasses="character"))

uniprot_length <- read.table("UniProt_len.txt",colClasses=c("character","numeric"))

PPI_edges <- read.table('PPI_edges.txt',colClasses="character",sep="\t")

cancer_types <- sapply(list.files('HQinterfaces_mutations/'),function(x) unlist(strsplit(x,split="_"))[1])

for(cancer in cancer_types){
    Mutations <- read.table(paste0('/data/Data/TCGA_Somatic_Mutations/Mutect2/TCGA.', cancer, '.mutect.NonSilent.maf'),
        skip=1,header=TRUE,colClasses="character",sep="\t")
    PPI_mutations <- read.table(paste0('HQinterfaces_mutations/', cancer, '_HQinterfaces_mutations.txt'),colClasses="character")
    uniprot_gene_mapping <- matrix(NA,nrow=nrow(PPI_edges),ncol=4)
    edge_pvalues <- matrix(NA,nrow=nrow(PPI_edges),ncol=6)
    for(i in 1:nrow(PPI_edges)){
        P1 <- PPI_edges[i,1]
        P1_res <- unlist(sapply(PPI_edges[i,4],function(x) strsplit(x,split=",")))
        P2 <- PPI_edges[i,2]
        P2_res <- unlist(sapply(PPI_edges[i,5],function(x) strsplit(x,split=",")))
        # P1
        mutation_number_1 <- sum(Mutations[,68]==P1)
        PPI_mutation_number_1 <- sum(PPI_mutations[,3]==P1 & PPI_mutations[,5] %in% P1_res)
        edge_pvalues[i,1:2] <- c(PPI_mutation_number_1,mutation_number_1)
        gene_len_1 <- uniprot_length[match(P1,uniprot_length[,1]),2]
        bind_len_1 <- length(P1_res)
        if(P1 %in% protein_coding_genes[,"uniprot_ids"]){
            uniprot_gene_mapping[i,1] <- protein_coding_genes[match(P1,protein_coding_genes[,"uniprot_ids"]),"symbol"]
            uniprot_gene_mapping[i,2] <- protein_coding_genes[match(P1,protein_coding_genes[,"uniprot_ids"]),"entrez_id"]
        }
        # P2
        mutation_number_2 <- sum(Mutations[,68]==P2)
        PPI_mutation_number_2 <- sum(PPI_mutations[,3]==P2 & PPI_mutations[,5] %in% P2_res)
        edge_pvalues[i,3:4] <- c(PPI_mutation_number_2,mutation_number_2)
        gene_len_2 <- uniprot_length[match(P2,uniprot_length[,1]),2]
        bind_len_2 <- length(P2_res)
        if(P2 %in% protein_coding_genes[,"uniprot_ids"]){
            uniprot_gene_mapping[i,3] <- protein_coding_genes[match(P2,protein_coding_genes[,"uniprot_ids"]),"symbol"]
            uniprot_gene_mapping[i,4] <- protein_coding_genes[match(P2,protein_coding_genes[,"uniprot_ids"]),"entrez_id"]
        }
        # Sig test
        p1 <- 1
        if(PPI_mutation_number_1 > 0 & mutation_number_1 > 0){
            p1 <- binom.test(x=PPI_mutation_number_1,n=mutation_number_1,p=bind_len_1/gene_len_1,alternative="greater")$p.value   
        }
        p2 <- 1
        if(PPI_mutation_number_2 > 0 & mutation_number_2 > 0){
            uniprot_gene_mapping[i,2] <- PPI_mutations[match(P2,PPI_mutations[,3]),2]
            p2 <- binom.test(x=PPI_mutation_number_2,n=mutation_number_2,p=bind_len_2/gene_len_2,alternative="greater")$p.value    
        }
        edge_pvalues[i,5] <- p1 * p2
    }
    result <- cbind(PPI_edges[,1:3],uniprot_gene_mapping,edge_pvalues)
    write.table(result,
                file=paste("./edge_pvalues/",cancer,"_edge_pvalues.txt",sep=""),
                row.names=FALSE,
                col.names=FALSE,
                sep="\t",
                quote=FALSE)
}

