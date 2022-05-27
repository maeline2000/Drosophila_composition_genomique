# Drosophila_composition_genomique
bio statistiques sur la composition génomique de 26 espèces du genre Drosophila
#TRAITEMENT DES EXONS
install.packages("BiocManager")
BiocManager::install("coRdon")

library("seqinr") # sequences manipulation
library("coRdon")
install.packages("dplyr")
library(tidyr)
library(dplyr)

#traitement des exons des espèces de Drosophiles

FAM_files_exons=Sys.glob("*.codingseq")
length(FAM_files_exons)

#fonction GC ############# seqinr ##################
GC_pct<-function(x){
  GC(getSequence(x))
}
#fonction GC3 
GC3_pct<-function(x){
  GC3(getSequence(x))
}

k=FAM_files_exons

#fonction CUB analysis ############# coRdon  ##################
CUB_analysis<-function(k){
  align_seqinr<-read.fasta(file=k, as.string = FALSE)
  align_coRdon<-readSet(file=k) 
  print(k)
  
  GC<-as.data.frame(do.call(rbind,lapply(align_seqinr, GC_pct)))
  GC3<-as.data.frame(do.call(rbind,lapply(align_seqinr, GC3_pct)))
  codon_table<-codonTable(align_coRdon)
  enc<- ENCprime(codon_table)
  MCB<-MCB(codon_table)
  scuo<-SCUO(codon_table)
  
  cub<- cbind(names(align_coRdon),enc,MCB,scuo,GC,GC3)
  names(cub)<-c("species","ENCprime","MCB","SCUO", "GC","GC3")
  
  #write.table(cub, file = paste(k,"_cub2.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
  return(cbind(cub,k))
}

cub_k<-lapply(FAM_files_exons[1:32], CUB_analysis) ## run the CUB analysis function on the files 

################ data analysis #######################
cub_table_sc<-do.call(rbind,cub_k[1:32]) #dataframe contening CUB stats for all genes 
names(cub_table_sc)<-c("species","ENCprime","MCB","SCUO", "GC","GC3","genes")
write.table(cub_table_sc,file = "cub_table_sc.txt",sep="\t",na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)

cub_table_sc<-read.table(file = "cub_table_sc.txt", dec = ".",h=TRUE)
head(cub_table_sc)
summary(cub_table_sc)
length(which(cub_table_sc$ENCprime<0))

###filtrer cub_table_sc :#filtre pour conserver les seqs orthologues présentes en une seule copie chez au moins la moitié des espèces étudiées
head(tab_filtre)
head(cub_table_sc)
o <- cub_table_sc
o1<-gsub(pattern = "g", replacement = "igg", x=o$species )
o1 <- as.data.frame(o1)
colnames(o1) <- c("o1")
o1<-separate(o1, o1, c('numéro_scaffold', 'gene'), sep = "ig")
gene_exon <- o1[,2]
gene_exon <- as.data.frame(gene_exon)
colnames(gene_exon) <- c("gene")
o2<-gsub(pattern = ".codingseq", replacement = "", x=o$genes )
o2 <- as.data.frame(o2)
colnames(o2) <- c("D_species_exon")
head(o2)

o3<-gsub(pattern = ".g", replacement = ".gg", x=o$species )
o3 <- as.data.frame(o3)
colnames(o3) <- c("gene_exon")
o3<-separate(o3, gene_exon, c('scaff','gene_exon'), sep = ".g")   #créer une colonne gene
gene_exon <- cbind(o2, o3)
gene_exon=mutate(gene_exon, gene=paste(D_species_exon, gene_exon, sep ="-") )
cub_table_sc <- cbind(cub_table_sc, gene_exon)
head(cub_table_sc)

#on regarde l'objet tab_filtre : obtenu avec le scrip dans le tri des introns
head(tab_filtre)  #on a la colonne gène en commun : jointure restreinte par "gene"
tab_group_sc <- inner_join(tab_filtre, cub_table_sc, by="gene")

head(tab_group_sc)
summary(tab_group_sc)


####représentation des taux de GC des séquences d'exons

###tableau avec les valeurs moyennes de CG par espece 
GC_sp_scf<-tab_group_sc
tab_GC_moy_par_espece_scf<-GC_sp_scf %>% group_by(D_species_exon) %>% summarise(GC= mean(GC))
head(tab_GC_moy_par_espece_scf)
plot_GCi_moy_par_espece_scf <- plot (tab_GC_moy_par_espece_scf, main = "GC moyen par espèce au niveau des séquences codantes \n chez différentes espèces de Drosophile  ",  cex.axis=0.75, las=2, xlab ="", ylab ="taux moyen de GC des séquences codantes")

####tableau avec les valeurs moyennes de CG3 par espece
GC3_sp_fsc<-tab_group_sc
tab_GC3_moy_par_espece_fsc<-GC3_sp_fsc %>% group_by(D_species_exon) %>% summarise(GC3= mean(GC3))
head(tab_GC3_moy_par_espece_fsc)
plot_GC3_moy_par_espece_fsc <- plot (tab_GC3_moy_par_espece_fsc, main = "GC3 moyen par espèce au niveau des séquences codantes \n chez différentes espèces de Drosophile  ",  cex.axis=0.75, las=2, xlab ="", ylab ="taux moyen de GC3 des séquences codantes")

###graphe du rapport GC / GC3 moyen par espèce
GC_SP_fsc<-tab_group_sc
TAB_GC_moy_par_espece_fsc<-GC_SP_fsc %>% group_by(D_species_exon) %>% summarise(GC3= mean(GC3), GC=mean(GC))
head(TAB_GC_moy_par_espece_fsc)
plot_TAB_GC_fsc <- plot (TAB_GC_moy_par_espece_fsc$D_species_exon, c(TAB_GC_moy_par_espece_fsc$GC / TAB_GC_moy_par_espece_fsc$GC3), main = "rapport entre le taux de GC moyen et le taux GC3 moyen \n par espèce au niveau des séquences codantes \n chez différentes espèces de Drosophile  ",  cex.axis=0.75, las=2, xlab ="", ylab ="taux moyen de GC/ taux moyen de GC3")

#pour les 3 espèces chez lesquels nous avons trouvé les taux de GC et de GC3 moyen les plus faibles, nous trouvons 
#un rapport GC/GC3 beaucoup plus proche de 1 que pour les espèces avec  des taux de GC plus élevés

###représentation du taux de GC en fonction des espèces
install.packages("ggplot2")
library(ggplot2)
GC_species_fsc <- tab_group_sc
GC_species_fsc$D_species_exon <- as.factor(GC_species_fsc$D_species_exon)
P_fsc <- ggplot(GC_species_fsc, aes(x=D_species_exon, y=GC, color=D_species_exon)) + labs(title = "Taux de GC moyen au niveau des séquences codantes en fonction des espèces de Drosophile", x="Espèce de Drosophile", y="Taux de GC moyen au niveau des séquences codantes")+ theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45)) + geom_boxplot() + stat_summary(fun = mean, geom="point", shape=20, size=4)
P_fsc

#TRAITEMENT DES INTRONS
###tableau avec les valeurs moyennes de CG par espece 
GCi_sp<-select(cub_table, species, GC)
GCi_sp

tab_GCi_moy_par_espece<-GCi_sp %>% group_by(species) %>% summarise(GC= mean(GC))
head(tab_GCi_moy_par_espece)
plot_GCi_moy_par_espece <- plot (tab_GCi_moy_par_espece, main = "GC moyen par espèce au niveau des introns")

###représentation du taux d'introns en fonction des espèces
library(ggplot2)
GCi_sp$species <- as.factor(GCi_sp$species)
p <- ggplot(GCi_sp, aes(x=species, y=GC, color=species)) + labs(title = "taux de GC des introns en fonction des espèces", x="espèce de Drosophile", y="taux de GC") + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45)) + stat_summary(fun = mean, geom="point", shape=20, size=4)
p



###GC moyen par groupe d'introns correspondant à un gène

GCi_genes<-select(cub_table, species, gene, GC, taille)
head(GCi_genes)
summary(GCi_genes)

GCi_genes$species <- as.factor(GCi_genes$species)
m <- ggplot(GCi_genes, aes(x=species, y=GC, color=species)) + labs(title = "Taux de GC des introns en fonction des espèces", x="espèce de Drosophile", y="taux de GC") + geom_boxplot() + stat_summary(fun = mean, geom="point", shape=20, size=4)
m 

GCi_genes_lm<-GCi_genes %>% group_by(species, gene) %>% summarise(GC= mean(GC), taille=sum(taille))
GCi_lm<-GCi_genes_lm %>% group_by(species, gene) %>% summarise(taille=sum(taille), GC=mean(GC))
head(GCi_lm)
summary(GCi_lm)
###pour affinis-g1t1 : GC moyen = 0,619 : on compare avec la valeur calculée manuellement pour vérifier : 
(0.4943 + 0.6875 + 0.6428 + 0.6518)/4 -> vm_affinis_t1g1         # = 0,6191

##considération de la taille des différentes séquences dans le calcul du GC/ (introns d'un gène)
GCi_lmp<-select(cub_table, species, gene, taille, GC)
GCi_lmp
GCi_lmp %>% group_by(gene) %>% summarise(weighted.mean(GC, taille))
GC_mp_théo <- ((0.4943*176 + 0.6875*112 + 0.6428*56 + 0.6518*135)/ (176+112+56+135)) #0,6012 = à la valeur obtenue pour g1-t1 dans GCi_lmp
#moyenne pondérée -> CG/gene différents : on refait les analyses avec ces nouvelles valeurs :
GCi_lmp$species <- as.factor(GCi_lmp$species)
mp <- ggplot(GCi_lmp, aes(x=species, y=GC, color=species)) + labs(title = "Taux de GC des introns en fonction des espèces", x="espèce de Drosophile", y="taux moyen de GC des introns par gène") + geom_boxplot() + stat_summary(fun = mean, geom="point", shape=20, size=4)
mp

GCi_gene_lmp_ff <- (GCi_lmp %>% group_by(gene, species) %>% summarise(weighted.mean(GC, taille)))
head(GCi_gene_lmp_filtre)
colnames(GCi_gene_lmp_filtre)[3] <- c("GC_lmp")
summary(GCi_gene_lmp_filtre)
mp_ggplot <- ggplot(GCi_gene_lmp_filtre, aes(x=species, y=GC_lmp, color=species)) + labs(title = "Taux de GC des introns filtrés en fonction des espèces", x="espèce de Drosophile", y="taux moyen de GC des introns par gène") + geom_boxplot() + stat_summary(fun = mean, geom="point", shape=20, size=4)
mp_ggplot

result_petits_introns <- by(TABLE_plot_petits_introns[,2:3], TABLE_plot_petits_introns$species, cor)
###### filter for outliers (ENC prime estimates <20)
filtered_cub_table<-subset(cub_table,ENCprime>20)
summary(filtered_cub_table)
write.table(filtered_cub_table,file = "filtered_cub_table.txt",sep="\t",na = "NA", dec = ".", row.names = FALSE,col.names = TRUE,append = FALSE)

#TRI DES INTRONS
#fichier de base
t2 <- read.csv("Dgp", header = FALSE)
head(t2)
t2[,1] #toute l'information est contenue dans cette colonne 1

t2 <- read.csv("Droso_gene_position", header = FALSE)

#définir l'objet comme un tableau et nommer un en tête de colonne "dgp"
t2 <- data.frame(t2)
colnames(t2) <- c("dgp")
head(t2)

#passer de 1 colonne à 2 : famille et identification du gène
install.packages(tidyr) 
library(tidyr)
t2<-separate(t2, dgp, c('famille','id_du_gene'), sep = " ")
head(t2)

famille<-t2[,1]
famille <- data.frame(famille)
colnames(famille) <- c("famille")
head(famille)

#créer une colonne espèce
t2<-separate(t2, id_du_gene, c('espece','idg', 'id', 'id2'), sep = "_")
head(t2)
#parfois un D devant le nom du genre drosophiles, d'autres fois non : on uniformise en supprimant les D, ou d correspondant à drosophila. : 
ti<-gsub(pattern = "dd", replacement = "popo", x=t2$espece )
ti<-gsub(pattern = "d", replacement = "", x=ti )
ti<-gsub(pattern = "popo", replacement = "D", x=ti )
ti<-gsub(pattern = "D", replacement = "", x=ti$species )
ti <- as.data.frame(ti)
colnames(ti) <- c("species")
head(ti)
summary (ti)

#création TABLE_INTERMEDIAIRE
#vérification que nos deux premières colonnes ont bien le même nombre de lignes
length(ti$species)
length(t2$famille)
#on les colle ensemble dans notre TABLE 
TABLE_2 <- (cbind(famille, ti))
head(TABLE_2)

#Sur t3 et t4 on extrait l'information du gène qu'on rajoute dans une 3ème colonne -> TABLE3 : 
t3 <-t2
t3 <- t2[,3:4]
t3=mutate(t3, id=paste(idg, id, sep ="_") )
t3 <- (t3[,2])
t3 <- as.data.frame(t3)
colnames(t3) <- c("ddd")

t3 <- gsub(pattern = "_NA", replacement = "", x=t3$ddd )
t3 <- as.data.frame(t3)
colnames(t3) <- c("ddd")

t3 <- gsub(pattern = "scaffold", replacement = "scaffold_", x=t3$ddd )
t3 <- as.data.frame(t3)
colnames(t3) <- c("ddd")

t3 <- gsub(pattern = "__", replacement = "_", x=t3$ddd )
t3 <- as.data.frame(t3)
colnames(t3) <- c("ddd")
head(t3)
length(t3$ddd)

t4 <- t3
t4 <- gsub(pattern = "scaffold_", replacement = "", x=t4$ddd )
t4 <- as.data.frame(t4)
colnames(t4) <- c("ddd")

#pour différencier le g de gène et celui de group qui décalent les colonnes -> Group
t4<-gsub(pattern = "grou", replacement = "Grou", x=t4$ddd )
t4 <- as.data.frame(t4)
colnames(t4) <- c("ddd")
head(t4)

t4<-gsub(pattern = ".g", replacement = ";;g", x=t4$ddd )
t4 <- as.data.frame(t4)
colnames(t4) <- c("ddd")

t4<-separate(t4, ddd, c('numéro_scaffold', 'gene', 'gene2'), sep = ";;")
head(t4)

length(t4$gene) #on obtient 122534 : même longueur que les autres colonnes de TABLE
g<- t4[,2]   #une colonne qui contient l'identification du gène
TABLE_3 <- (cbind(famille, ti, g))
head(TABLE_3)
length(TABLE_3[,3])

###on extrait l'info : scaffold : on la met au même format que dans cub.table : scaffold-x
t5 <- t4
t5=mutate(t5, scaffold=paste(gene2, numéro_scaffold, sep ="-") )
t5 <- (t5[,4])
t5 <- as.data.frame(t5)
colnames(t5) <- c("scaffold")
head(t5)

t5 <- gsub(pattern = "NA", replacement = "scaffold", x=t5$scaffold )
t5 <- as.data.frame(t5)
colnames(t5) <- c("scaffold")
head(t5)

length(t5$scaffold) #on obtient 122534 : même longueur que les autres colonnes de TABLE
scaffold<- t5[,1]   #une colonne qui contient l'identification du gène
TABLE_4 <- (cbind(famille, ti, scaffold, g))
head(TABLE_4)

###on créer une colonne gène identique à celle de l'autre fichier : espece-gx.tx
t6 <- TABLE_3
t6=mutate(t6, gene=paste(species, g, sep ="-") )
t6 <- (t6[,4])
t6 <- as.data.frame(t6)
colnames(t6) <- c("gene")
head(t6)
length(t6[,1])

length(t6$gene) #toujours la même longueure de colonne -> ok
gene  <- t6$gene
gene  <- as.data.frame(gene)
TABLE_5 <- (cbind(famille, ti, scaffold, gene))
head(TABLE_5)

#maintenant que nous avons les colonnes qui nous intéressent dans nos deux tableaux de données 
#on filtre les gènes des espèces que nous gardons : 
tab_filtre <- (cbind(famille, scaffold, gene))
head(tab_filtre)
head(cub_table)

#on fait une jointure restreinte pour ne garder que les les lignes correspondantes entre les deux tableaux 
#celles présentes dans un seul des deux data.frame ne sont pas conservées
tab_group <- inner_join(tab_filtre, cub_table, by="gene")
head(tab_group)
summary(tab_group)
t<-tab_group %>% group_by(gene) %>% summarise(GC= mean(GC))
head(t)
summary(t)##NB : aucune séquences d'algonquin n'a été conservée : le jeu de données à partir duquel nous avons filtré était-il complet?

#le data.frame tab_group conserve uniquement les séquences que nous voulons conserver pour l'analyse : représentation du jjd conservé : 
### représentation_seq_introns_filtrés : on rajoute _f à la suite du nom des fichiers ou les séquences ont été filtré
#tableau avec les valeurs moyennes de CG par espece 
GCi_sp_f<-select(tab_group, species, GC)
head(GCi_sp_f)
summary(GCi_sp_f)
tab_GCi_moy_par_espece_f<-GCi_sp_f %>% group_by(species) %>% summarise(GC= mean(GC))
head(tab_GCi_moy_par_espece_f)
#on peut comparer avec la valeur obtenue avant d'effetuer un filtre pour affinis : 
head(tab_GCi_moy_par_espece)
#avant le filtre : GC moyen =0.421 pour affinis, et après le filtre : taux de GC moyen = 0.414

plot_GCi_moy_par_espece_f <- plot (tab_GCi_moy_par_espece_f, main = "GC moyen par espèce au niveau des introns filtrés" )

###représentation du taux d'introns en fonction des espèces
library(ggplot2)
GCi_sp_f$species <- as.factor(GCi_sp_f$species)
p_f <- ggplot(GCi_sp_f, aes(x=species, y=GC, color=species)) + labs(title = "taux de GC des introns filtrés en fonction des espèces", x="espèce de Drosophile", y="taux de GC") + geom_boxplot() + stat_summary(fun.y = mean, geom="point", shape=20, size=4)
p_f

###GC moyen par groupe d'introns filtrés correspondant à un gène : on considère les moyennes pondérées des longueurs des séquences :
library(tidyr)
library(dplyr)
#considération de la taille des différentes séquences dans le calcul du GC/ (introns d'un gène)
GCi_gene_lmp_f<-select(tab_group, species, gene, taille, GC)
summary(GCi_gene_lmp_f)

GCi_gene_lmp_filtre <- (GCi_gene_lmp_f %>% group_by(gene, species) %>% summarise(weighted.mean(GC, taille)))
head(GCi_gene_lmp_filtre)
colnames(GCi_gene_lmp_filtre)[3] <- c("GC_lmp")
summary(GCi_gene_lmp_filtre)
mp <- ggplot(GCi_gene_lmp_filtre, aes(x=species, y=GC_lmp, color=species)) + labs(title = "Taux de GC des introns filtrés en fonction des espèces", x="espèce de Drosophile", y="taux moyen de GC des introns par gène") + geom_boxplot() + stat_summary(fun.y = mean, geom="point", shape=20, size=4)
mp

#TABLE_GCi_GC3 PAR ESPÈCE
###ici on a les GCi par gène et par espèce : et si on essaie de rajouter une colonne GC3 par gène pour les différentes espèces? 
#on a : cub_table_sc_f : qui contient les gènes filtrés ainsi que leur taux de GC3 : 
head(cub_table_sc)
gene_TABLE_plot <- cub_table_sc[,11]
gene_TABLE_plot <- data.frame(gene_TABLE_plot)
colnames(gene_TABLE_plot) <- c("gene")

GC3_TABLE_plot <- cub_table_sc[,6]
GC3_TABLE_plot <- data.frame(GC3_TABLE_plot)
colnames(GC3_TABLE_plot) <- c("GC3")

species_TABLE_plot <- cub_table_sc[,8]
species_TABLE_plot <- data.frame(species_TABLE_plot)
colnames(species_TABLE_plot) <- c("species")

TABLE_plot_A <- (cbind(species_TABLE_plot, gene_TABLE_plot, GC3_TABLE_plot))

head(GCi_gene_lmp_filtre)
gene_TABLE_plot_B <- GCi_gene_lmp_filtre[,1]
gene_TABLE_plot_B <- data.frame(gene_TABLE_plot_B)
colnames(gene_TABLE_plot_B) <- c("gene")

GCi_TABLE_plot <- GCi_gene_lmp_filtre[,3]
GCi_TABLE_plot <- data.frame(GCi_TABLE_plot)
colnames(GCi_TABLE_plot) <- c("GCi")

TABLE_plot_B <- (cbind(gene_TABLE_plot_B, GCi_TABLE_plot))

TABLE_plot <- inner_join(TABLE_plot_A, TABLE_plot_B, by="gene")
head(TABLE_plot)

#à partir du tableur : TABLE_plot on essaie de tracer les corrélations GC3 _ GCi par espèce : fonction ggplot2?
install.packages("tidyverse")
library(tidyverse)
plot_correlation_group <- ggplot(data = TABLE_plot) + facet_wrap(vars(species)) + geom_point(aes(x = GC3, y = GCi, color = species)) + labs(title = "taux de GCi en fonction du taux de GC3 par espèce de Drosophile")
plot_correlation_group
result <- by(TABLE_plot[,3:4], TABLE_plot$species, cor)
result

TABLE_correlation <- TABLE_plot %>% group_by(species) %>% summarize(cor.test=cor.test(GCi, GC3))
TABLE_correlation

#CORRELATION POUR LES PETITS INTRONS (MOINS DE 80 PB)
#on rajoute une colonne au tableau correspondant à la longueur des séquences
head(cub_table)
petits_introns <- cub_table
petits_introns=filter(petits_introns, length<80)
head(petits_introns)  #table contenant uniquement les valeurs de GC des petits introns
petits_introns <- petits_introns[,2:6]
tab_petits_introns <- inner_join(tab_filtre, petits_introns, by="gene") #on garde uniquement les introns appartenant aux orthologues
tab_petits_introns <- inner_join(tab_petits_introns, TABLE_plot_A, by="gene") #on rajoute une colonne correspondant au GC3 des gènes
head(tab_petits_introns)
i <- tab_petits_introns[,8]
i <- as.data.frame(i)
colnames(i) <- c("species")
ii <- tab_petits_introns[,5]
ii <- as.data.frame(ii)
colnames(ii) <- c("GCi")
iii <- tab_petits_introns[,9]
iii <- as.data.frame(iii)
colnames(iii) <- c("GC3")
TABLE_plot_petits_introns <- cbind(i,ii,iii)

plot_correlation_group_petits_introns <- ggplot(data = TABLE_plot_petits_introns) + facet_wrap(vars(species)) + geom_point(aes(x = GC3, y = GCi, color = species)) + labs(title = "taux de GC au niveau des introns de taille inférieur à 80 pb \n en fonction du taux de GC3 par espèce de Drosophile ")
plot_correlation_group_petits_introns

###petits introns
head(petits_introns)
pi<-gsub(pattern = "-g", replacement = "-ggg", x=petits_introns$id )
pi <- as.data.frame(pi)
colnames(pi) <- c("id")
pi<-separate(pi, id, c('idd','gene'), sep = "-gg")
petits_introns <- cbind(petits_introns, pi)
head(petits_introns)
TAB_pi<-petits_introns %>% group_by(gene, species) %>% summarise(GCi= mean(GCi))
head(TAB_pi)
pii=mutate(TAB_pi, gene=paste(species, gene, sep ="-") )
pii <- as.data.frame(pii)
pii <- cbind(pii$gene, pii$GCi)
pii <- as.data.frame(pii)
colnames(pii) <- c("gene","GCi")
tab_data <- inner_join(pii, cub_table_sc, by="gene")
head(tab_data)

head(TABLE_plot_petits_introns)
result_petits_introns <- by(TABLE_plot_petits_introns[,2:3], TABLE_plot_petits_introns$species, cor)
result_petits_introns

TABLE_correlation_petits_introns <- TABLE_plot_petits_introns %>% group_by(species) %>% summarize(cor.test=cor.test(GCi, GC3))
TABLE_correlation_petits_introns
affinis_seq <- filter(TABLE_plot_petits_introns, species=="affinis")
ambigua_seq <- filter(TABLE_plot_petits_introns, species=="ambigua")
athabasca_seq <- filter(TABLE_plot_petits_introns, species=="athabasca")
bifasciata_seq <- filter(TABLE_plot_petits_introns, species=="bifasciata")
ercepae_seq <- filter(TABLE_plot_petits_introns, species=="ercepae")
helvetica_seq <- filter(TABLE_plot_petits_introns, species=="helvetica")
lucipennis_seq <- filter(TABLE_plot_petits_introns, species=="lucipennis")
malerkotliana_seq <- filter(TABLE_plot_petits_introns, species=="malerkotliana")
guanche_seq <- filter(ttg, species=="guanche")
imaii_seq <- filter(ttg, species=="imaii")
me07rcatorum_seq <- filter(ttg, species=="me07rcatorum")
microlabis_seq <- filter(ttg, species=="microlabis")
nebulosa_seq <- filter(ttg, species=="nebulosa")
phaeopleura_seq <- filter(ttg, species=="phaeopleura")
psaltans_seq <- filter(ttg, species=="nebulosa")
pstakahashi_seq <- filter(ttg, species=="pstakahashi")

head(tab_data)
as.numeric(tab_data$GCi)
r <- as.data.frame.numeric(tab_data$GCi)
colnames(r) <- c("GCi")
rr <- as.data.frame.numeric(tg$GC3)
colnames(rr) <- c("GC3")
rrr <- as.data.frame(tg$D_species_exon)
colnames(rrr) <- c("species")
td <- cbind(rrr, r, rr)
head(td)
ttg <- transform(td, GCi = as.numeric(GCi))

cor.test(guanche_seq$GC3, guanche_seq$GCi)
cor.test(imaii_seq$GC3, imaii_seq$GCi)
cor.test(me07rcatorum_seq$GC3, me07rcatorum_seq$GCi)
cor.test(microlabis_seq$GC3, microlabis_seq$GCi)
cor.test(nebulosa_seq$GC3, nebulosa_seq$GCi)
cor.test(phaeopleura_seq$GC3, phaeopleura_seq$GCi)
cor.test(psaltans_seq$GC3, psaltans_seq$GCi)
cor.test(pstakahashi_seq$GC3, pstakahashi_seq$GCi)

#ARBRE PHYLOGÉNÉTIQUE
install.packages("ape")
library("ape")
tree <- read.tree(file="root_droso_MF_for_publication")
tree1 <- plot(tree, main="Arbre phylogenetique de 26 espèces du genre Drosophila", cex=0.6, edge.width=1.5, edge.color = "deepskyblue4")

working_tree <- tree
class(working_tree)
working_tree <- as.phylo(tree)
tip <- c("D.varians*", "D.pseudoananassae*", "D.bipectinata", "D.paulistorum*", "D.pallidosa*", "D.ananassae", "D.kikkawai", "D.serrata", "D.diplacantha*", "D.ficusphila", "D.elegans", "D.rhopaloa", "D.biarmipes", "D.suzukii", "D.pseudotakahashi*", "D.takahashii", "D.eugracilis", "D.erecta", "D.yakuba", "D.teissieri*", "D.melanogaster", "D.sechellia", "D.simulans", "D.mauritiana", "D.madeirensis*", "D.subobscura", "D.subobscura*", "D.tsukubaensis*","D.lowei", "D.miranda", "D.persimilis","D.pseudoobscura", "D.willistoni", "D.busckii", "D.nasuta", "D.albomicans", "Z.indianus*",  "Z.gabonicus*", "Z.africanus*", "D.grimshawi","D.lacertosa","D.nigromelanica", "D.micromelanica", "D.melanica", "D.robusta","D.montana", "D.virilis","D.novamexicana","D.americana", "D.hydei", "D.navojoa", "D.mojavensis","D.arizonae")
tree <- plot(drop.tip(working_tree, tip),  main="Arbre phylogenetique de 26 espèces du genre Drosophila", cex=0.6, edge.width=1.5, edge.color = "deepskyblue4")
