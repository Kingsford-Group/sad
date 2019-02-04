library(polyester)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

fasta = readDNAStringSet(args[1])
readlen = 100
readspertx = round(1 * width(fasta) / readlen)
fold_changes = matrix(rep(1, 2*length(fasta)), nrow=length(fasta))

theoretical <- read.table(args[2], header=FALSE, sep="\t")
colnames(theoretical) <- c("Name", "NumReads")

theoretical <- theoretical[match(names(fasta), theoretical$Name), ]
for (i in 1:length(fasta)) {
	readspertx[i] = max(round(theoretical[i, "NumReads"]), 1)*as.numeric(args[4])
}

simulate_experiment(args[1], reads_per_transcript=readspertx, num_reps=c(1,1), fold_changes=fold_changes, error_model='illumina5', 
	bias='cdnaf', outdir=args[3])


# t1<-read.table('rank_singleratio.txt', header=FALSE, sep="\t")
# t2<-read.table('tmpGM12878_emd.txt', header=FALSE, sep="\t")
# colnames(t1)<-c("Name", "ratio")
# colnames(t2)<-c("Name", "EMD")
# overlap<-intersect(t1$Name, t2$Name)
# t1<-t1[match(overlap, t1$Name),]
# t2<-t2[match(overlap, t2$Name),]
# df<-data.frame(t1$Name, t1$ratio, t2$EMD)
# colnames(df)<-c("Name", "ratio", "EMD")
# corr<-cor(df$ratio, df$EMD)
# ggplot(df)+geom_point(aes(x=EMD, y=ratio))+geom_text(aes(x=0.4, y=0.9, label=paste("pearson cor = ", round(corr, digits=3), sep="")))