### TOOLS
library(ape)
library(mgcv)
binarize <- function(x, maxn) as.numeric(seq(maxn) %in% x)  
  

### GENE TREES & BIPARTITION PROBABILITIES
# complete list of species
species <- c("Aethomys_hindei", "Aethomys_chrysophilus", "Aethomys_kaiseri", "Aethomys_nyikae", 
	"Arvicanthis_blicki", "Arvicanthis_nairobae", "Arvicanthis_niloticus", "Arvicanthis_rufinus", "Arvicanthis_somalicus", 
	"Dasymys_cf_incomtus", "Dephomys_defua", "Desmomys_harringtoni", "Golunda_ellioti", 
	"Grammomys_dolichurus_group", "Grammomys_macmillani_group", "Grammomys_poensis_group", "Grammomys_selousi_group", "Grammomys_surdaster_group", 
	"Hybomys_cf_lunarisnsp", "Hybomys_planifrons", "Hybomys_trivirgatus", "Hybomys_univittatus", "Lamottemys_okuensis", 
	"Lemniscomys_barbarus", "Lemniscomys_bellieri", "Lemniscomys_rosalia", "Lemniscomys_striatus", "Micaelamys_namaquensis", 
	"Millardia_meltada", "Mylomys_dybowskii", "Oenomys_hypoxanthus", "Otomys_sp_2", "Otomys_typus_sl", "Pelomys_fallax", "Praomys_delectorum", 	"Rhabdomys_dilectus", "Stenocephalemys_albocaudata", "Stochomys_longicaudatus", "Thallomys_paedulcus", "Thamnomys_kempi")

# burnin
burnin <- 501

# path to a folder with MrBayes outputs for individual loci
# MrBayes outputs = tree files ('.t') & log files ('.p')
path <- "mbresults_arvicanthini"
# listing the files and splitting them according to loci
# the files are named LOCUSNAME.run#.t & LOCUSNAME.run#.p
files <- list.files(path, full.names=TRUE)
loci <- gsub("\\.run[[:digit:]]\\.[[:alpha:]]$", "", gsub(paste(path, "/", sep=""), "", files))
loci <- setNames(split(files, loci), unique(loci))

# variable to store posterior probabilities of sampled bipartitions
PP <- setNames(vector("list", length(loci)), names(loci))

# loop calculating PP for every locus
for (i in seq_along(loci)) {
	treefiles <- loci[[i]][grepl("\\.t", loci[[i]])]
	trees <- lapply(lapply(treefiles, ape::read.nexus), "[", -seq(burnin))
	tips <- trees[[1]][[1]]$tip.label
	part <- lapply(trees, ape::prop.part, check.labels=FALSE)
	counts <- unlist(lapply(part, attr, "number"))
	part <- unlist(part, recursive=FALSE)
	part <- lapply(lapply(part, function(x) tips[x]), match, table=species)
	part <- apply(sapply(part, binarize, maxn=length(species)), 2, paste, collapse="")
	PP[[i]] <- sapply(split(counts, part), sum) / sum(sapply(trees, length))
}
saveRDS(PP, "PP_190905.rds")
PP <- readRDS("PP_190905.rds")
PP <- readRDS("data/PP_180808")

### SELECTION OF LOCI
# filtering out loci with incomplete set of species 
PP <- PP[sapply(PP, function(x) any(grepl("1", names(x)) & !grepl("0", names(x))))]
bipartitions <- sort(unique(unlist(lapply(PP, names))))

# calculation of the strength & diversity of phylogenetic signal
L <- matrix(0, length(PP), length(bipartitions), dimnames=list(names(PP), bipartitions))
for (i in seq_along(PP)) {
	L[i, names(PP[[i]])] <- PP[[i]]
}
P <- L %*% t(L) / length(species)
ord <- order(diag(P), decreasing=TRUE)
P <- P[ord, ord]
k <- seq(nrow(P))[-1]
S <- sapply(k, function(i) mean(diag(P)[seq(i)]))
D <- 1 - sapply(k, function(i) mean(P[seq(i),seq(i)][upper.tri(diag(i))]))
Q <- S * D

# smoothing of Q over the increasing size of the corresponding subset of loci 
smooth <- predict(mgcv::gam(Q ~ s(k, bs="cr")), type="response")
# calculation of threshold
thr <- which.max(diff(smooth, differences=3)) + 3 + 1
# plotting of the result
plot(k, Q, pch=16, col=3)
lines(k, smooth, lwd=2)
points(thr - 1, smooth[thr - 1], pch=16, col=2)
# export of the names of selected loci
loci <- rownames(P)[seq(thr)]
writeLines(loci, "selected_loci.txt")

