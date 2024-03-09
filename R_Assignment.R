library(tidyverse)
library(dplyr)
library(RColorBrewer)

#Reading in data
fang <- read.delim("fang_et_al_genotypes.txt")
snp <- read.delim("snp_position.txt")


#subsetting Maize data
maize <- filter(fang, Group=="ZMMIL" | Group=="ZMMLR" | Group=="ZMMMR")
write.csv(maize, "maize.csv")
#transposing
trans_maize <- data.frame(t(maize))
#Adding column name
trans_maize <- cbind(rownames(trans_maize), trans_maize)
rownames(trans_maize) <- NULL
colnames(trans_maize) <- c(names(trans_maize))
colnames(trans_maize)[1] <- "Sample_ID"
#Checking
names(trans_maize)
write.csv(trans_maize, "trans_maize.csv")

#merging with SNP
maize_merged <- merge(snp, trans_maize,by.x="SNP_ID", by.y = "Sample_ID")
write.csv(maize_merged, "merged_maize.csv")
#pulling chromsome files, sorting
for (i in 1:10) {
  subset <- filter(maize_merged, Chromosome==i)
  subset2 <- subset[order(as.numeric(subset$Position)),]
  write.csv(subset2, paste("maize", i, "csv", sep = "."))
}

#replacing question marks
maize_merged[-1] <- lapply(maize_merged[-1],gsub,pattern ="\\?",replacement ="-")
write.csv(maize_merged, "maize_merged.csv")
#ulling chromosome files, sorting
for (i in 1:10) {
  subset <- filter(maize_merged, Chromosome==i)
  subset2 <- subset[order(as.numeric(subset$Position), decreasing = TRUE),]
  write.csv(subset2, paste("sub_maize", i, "csv", sep = "."))
}


#subsetting Teosinte data
teosinte <- filter(fang, Group=="ZMPBA" | Group=="ZMPIL" | Group=="ZMPJA")

#transposing
trans_teosinte <- data.frame(t(teosinte))
#Adding column name
trans_teosinte <- cbind(rownames(trans_teosinte), trans_teosinte)
rownames(trans_teosinte) <- NULL
colnames(trans_teosinte) <- c(names(trans_teosinte))
colnames(trans_teosinte)[1] <- "Sample_ID"

names(trans_teosinte)

#merging with SNP
teosinte_merged <- merge(snp, trans_teosinte,by.x="SNP_ID", by.y = "Sample_ID")


#pulling chromsome files, sorting
for (i in 1:10) {
  subset <- filter(teosinte_merged, Chromosome==i)
  subset2 <- subset[order(as.numeric(subset$Position)),]
  write.csv(subset2, paste("teosinte", i, "csv", sep = "."))
}

#replacing question marks
teosinte_merged[-1] <- lapply(teosinte_merged[-1],gsub,pattern ="\\?",replacement ="-")

#pulling chromosome files, sorting
for (i in 1:10) {
  subset <- filter(teosinte_merged, Chromosome==i)
  subset2 <- subset[order(as.numeric(subset$Position), decreasing = TRUE),]
  write.csv(subset2, paste("sub_teosinte", i, "csv", sep = "."))
}
trans_fang <- data.frame(t(fang))
trans_fang <- cbind(rownames(trans_fang), trans_fang)
rownames(trans_fang) <- NULL
colnames(trans_fang) <- c(names(trans_fang))
colnames(trans_fang)[1] <- "Sample_ID"
all_merged <- merge(snp, trans_fang, by.x="SNP_ID", by.y = "Sample_ID")

write.csv(all_merged,"trans_fang.csv")
#subsetting site columns
merged_maize_x <- maize_merged[ , grepl("X",names(maize_merged))]
#adding new column for homozygous or not
homozygous <- c("A/A", "C/C", "T/T", "G/G")
merged_maize_x <- merged_maize_x %>%
 mutate(Homozygous = pmap(.,~any(homozygous %in% .)))

merged_maize_x <- apply(merged_maize_x,2,as.character)
write.csv(merged_maize_x, "maize_x.csv")
#plot with counts
ggplot(merged_maize_x, aes(x = Homozygous)) + geom_bar()+ geom_text(stat='count',aes(label=after_stat(count)), vjust=1)

#subsetting site columns
merged_teo_x <- teosinte_merged[ , grepl("X",names(teosinte_merged))]
#adding new column for homozygous or not
merged_teo_x <- merged_teo_x %>%
  mutate(Homozygous = pmap(.,~any(homozygous %in% .)))



merged_teo_x <- apply(merged_teo_x,2,as.character)
write.csv(merged_teo_x, "teo_x.csv")
#plot with counts
ggplot(merged_teo_x, aes(x = Homozygous)) + geom_bar() + geom_text(stat='count',aes(label=after_stat(count)), vjust=1)


#New plot
merged_maize_x <- as.data.frame(merged_maize_x)
maize_x <- data.frame(merged_maize_x$Homozygous)

merged_teo_x <- as.data.frame(merged_teo_x)
teo_x <- data.frame(merged_teo_x$Homozygous)

x_data <- data.frame(maize_x,teo_x)
colnames(x_data) <- c('maize', 'teosinte')
x_data <- as.data.frame(x_data)
x_data_2 <- pivot_longer(x_data,cols=c('maize', 'teosinte'),names_to = 'group', values_to = 'Homozygous')

ggplot(x_data_2,aes(fill=Homozygous,x=group)) + geom_bar() + geom_text(stat='count',aes(label=after_stat(count)),size = 3, position = position_stack(vjust = 0.5))



maize_merged_chrom <- data.frame(maize_merged$Chromosome)
teosinte_merged_chrom <- data.frame(teosinte_merged$Chromosome)

chrom_data <- data.frame(maize_merged_chrom, teosinte_merged_chrom)
colnames(chrom_data) <- c('maize', 'teosinte')
chrom_data_2 <- pivot_longer(chrom_data, cols=c('maize', 'teosinte'), names_to = 'group', values_to = 'Chromosome')

ggplot(chrom_data_2, aes(fill=Chromosome,x = group)) + geom_bar(position = "dodge")


snp_pos <- data.frame(snp$Chromosome,snp$Position)
colnames(snp_pos) <- c('Chromosome', 'Position')
snp_pos$Position <- as.numeric(snp_pos$Position)
snp_pos <- snp_pos %>%
  mutate(pos_binned = cut(Position,5));
snp_pos <- snp_pos[!snp_pos$Chromosome %in% c('unknown','multiple'),]
snp_pos <- na.omit(snp_pos)


pos_plot <- ggplot(data=snp_pos) + geom_bar(mapping=aes(x=pos_binned, fill = Chromosome))
pos_plot + theme(axis.text.x= element_text(size = 3)) + xlab("Position") + scale_fill_brewer(palette = 'Spectral')
