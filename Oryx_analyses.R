# load libraries
library(radiator)
library(dartR)
library(readr)
library(stringr)
library(gsubfn)
library(ggplot2)
library(gplots)
library(plotly)
library(hierfstat)
library(poppr)
library(pegas)
library(sequoia)
library(kinship2)
library(bigsnpr)
library(plinkQC)

########### INPUT AND FILTERING ###########
# Individuals that were sampled twice that will be removed from the final database
#          NAs       ID      NAs_percent
# OXL487b  2244  OXL487b          57
# OXL581b  1313  OXL581b          33
# OXL505b  1169  OXL505b          30
# OXL505a   631  OXL505a          16
# individual 471 is the same individual as 472 and was mislabelled. so it is exluded from analyses
individuals_to_remove <- c("OXL487b","OXL581b","OXL505b","OXL505a","OXL471")

oryx_vcf_radiator <- radiator::read_vcf(
  data = "2_5.vcf",
  strata = "strata_oryx.txt",
  filter.individuals.missing = 0.5,
  filter.common.markers = TRUE,
  filter.mac = 3,
  filter.genotyping = 0.1,
  filter.short.ld = "mac",
  filter.long.ld = NULL,
  verbose = TRUE
)

oryx_gl <- write_genlight(oryx_vcf_radiator,dartr = T)

oryx_gl <- gl.compliance.check(oryx_gl)
oryx_gl <- gl.recalc.metrics(oryx_gl)
ploidy(oryx_gl) <- 2

oryx_gl <- gl.drop.ind(oryx_gl, ind.list=individuals_to_remove)
oryx_gl <- gl.filter.hwe(oryx_gl)
# Removing "OXL" from names 
oryx_gl$ind.names <- substring(oryx_gl$ind.names,4)
oryx_gl$ind.names <- substring(oryx_gl$ind.names,1,3)
# separating the words in the names of populations 
oryx_gl$pop <-as.factor(sub("_"," ",oryx_gl$pop))
# individual 454 was mislabelled and was changed from WWR_Oman to WWR-Mix
ind_454 <- which(oryx_gl$ind.names=="454")
oryx_gl$pop[ind_454] <- as.factor("WWR Mix")
oryx_gl <- gl.recalc.metrics(oryx_gl)

# subsetting populations
samples_oryx <- read.csv("oryx_samples_radiator.csv")
oryx_gl$other$info <- samples_oryx
levels(pop(oryx_gl))
Mix <- gl.keep.pop(oryx_gl, pop.list="WWR Mix",mono.rm = F ) 
Oman <- gl.keep.pop(oryx_gl, pop.list="WWR Oman",mono.rm = F ) 
UAE <- gl.keep.pop(oryx_gl, pop.list="WWR UAE",mono.rm = F) 
historical <- gl.keep.pop(oryx_gl, pop.list=c("Oman founder","USA founder"),mono.rm = F ) 
without_historical <- gl.keep.pop(oryx_gl, pop.list=c("WWR Mix","WWR Oman","WWR UAE"),mono.rm = F ) 
Oman_UAE <-  gl.keep.pop(oryx_gl, pop.list=c("WWR Oman","WWR UAE"),mono.rm = F ) 
##################################
######### STRUCTURE ##############
##################################
######### ALL ####################
##################################
write.csv(cbind(oryx_gl$ind.names,as.character(oryx_gl$pop)),file = "order_structure.csv")
order_structure <- read_csv("order_structure_ordered.csv")
oryx_structure <- oryx_gl
oryx_structure$other$info <- order_structure
matrix_structure <- as.data.frame(oryx_structure)
matrix_structure[,(ncol(matrix_structure)+1)] <- oryx_structure$other$info$pop_structure
columns <- ncol(matrix_structure)
# sorting individuals based on previous structure analyses that were aimed in sorting the individuals by q in distruct
matrix_structure <- matrix_structure[order(matrix_structure[,columns]),]
matrix_structure_2 <- as.matrix(matrix_structure)
matrix_structure_2 <- matrix_structure_2[,-columns]
matrix_structure_2 <- as.genlight(matrix_structure_2)
ploidy(matrix_structure_2) <- 2
gl2faststructure(matrix_structure_2, outfile = "oryx_all_2.str",outpath=getwd())
pop_names_sorted_structure <- as.data.frame(cbind(oryx_structure$other$info$pop_structure,oryx_structure$ind.names,as.character(oryx_structure$pop)))
# this is the populations names to plot in distruct
colnames(pop_names_sorted_structure) <- c("pop_structure","ID","pop")
pop_names_sorted_structure$pop_structure <- as.numeric(as.character(pop_names_sorted_structure$pop_structure))
pop_names_sorted_structure <- pop_names_sorted_structure[order(pop_names_sorted_structure$pop_structure),]
write.table(as.data.frame(pop_names_sorted_structure$ID),file = "pops_per_individual_all.txt",quote = F,col.names = F,row.names = F)
##################################
###### WITHOUT HISTORICAL ########
##################################
matrix_structure_wo_hist <- as.data.frame(without_historical)
matrix_structure_wo_hist[,(ncol(matrix_structure_wo_hist)+1)] <- without_historical$other$info$pop_structure
columns_wo_hist <- ncol(matrix_structure_wo_hist)
# sorting individuals based on previous structure analyses that were aimed in sorting the individuals by q in distruct
matrix_structure_wo_hist <- matrix_structure_wo_hist[order(matrix_structure_wo_hist[,columns_wo_hist]),]
matrix_structure_2_wo_hist <- as.matrix(matrix_structure_wo_hist)
matrix_structure_2_wo_hist <- matrix_structure_2_wo_hist[,-columns_wo_hist]
matrix_structure_2_wo_hist <- as.genlight(matrix_structure_2_wo_hist)
ploidy(matrix_structure_2_wo_hist) <- 2
gl2faststructure(matrix_structure_2_wo_hist, outfile = "without_historical.str", probar = TRUE,outpath=getwd())
pop_names_sorted_structure_wo_hist <- as.data.frame(cbind(without_historical$other$info$pop_name_structure,without_historical$other$info$pop_structure,without_historical$ind.names))
# this is the populations names to plot in distruct
colnames(pop_names_sorted_structure_wo_hist) <- c("pop_structure","order_structure","ID")
pop_names_sorted_structure_wo_hist$order_structure <- as.numeric(as.character(pop_names_sorted_structure_wo_hist$order_structure))
pop_names_sorted_structure_wo_hist <- pop_names_sorted_structure_wo_hist[order(pop_names_sorted_structure_wo_hist$order_structure),]
# the following names should be pasted first into excel so there is not mistake with cumplak
cat(as.character(pop_names_sorted_structure_wo_hist$pop_structure), sep='\n')
cat(as.character(pop_names_sorted_structure_wo_hist$ID), sep='\n')
##################################
###### Oman_UAE #################
##################################
matrix_structure_Oman_UAE <- as.data.frame(Oman_UAE)
matrix_structure_Oman_UAE[,(ncol(matrix_structure_Oman_UAE)+1)] <- Oman_UAE$other$info$pop_structure
columns_Oman_UAE <- ncol(matrix_structure_Oman_UAE)
# sorting individuals based on previous structure analyses that were aimed in sorting the individuals by q in distruct
matrix_structure_Oman_UAE <- matrix_structure_Oman_UAE[order(matrix_structure_Oman_UAE[,columns_Oman_UAE]),]
matrix_structure_2_Oman_UAE <- as.matrix(matrix_structure_Oman_UAE)
matrix_structure_2_Oman_UAE <- matrix_structure_2_Oman_UAE[,-columns_Oman_UAE]
matrix_structure_2_Oman_UAE <- as.genlight(matrix_structure_2_Oman_UAE)
ploidy(matrix_structure_2_Oman_UAE) <- 2
gl2faststructure(matrix_structure_2_Oman_UAE, outfile = "Oman_UAE.str", probar = TRUE,outpath=getwd())
pop_names_sorted_structure_Oman_UAE <- as.data.frame(cbind(Oman_UAE$other$info$pop_name_structure,Oman_UAE$other$info$pop_structure,Oman_UAE$ind.names))
# this is the populations names to plot in distruct
colnames(pop_names_sorted_structure_Oman_UAE) <- c("pop_structure","order_structure","ID")
pop_names_sorted_structure_Oman_UAE$order_structure <- as.numeric(as.character(pop_names_sorted_structure_Oman_UAE$order_structure))
pop_names_sorted_structure_Oman_UAE <- pop_names_sorted_structure_Oman_UAE[order(pop_names_sorted_structure_Oman_UAE$order_structure),]
# the following names should be pasted first into excel so there is not mistake with cumplak
cat(as.character(pop_names_sorted_structure_Oman_UAE$pop_structure), sep='\n')
cat(as.character(pop_names_sorted_structure_Oman_UAE$ID), sep='\n')

##################################
###### ESTIMATING K #################
##################################
files_structure <- list.files(path = paste0(getwd(),"/clumpak_oryx_all_ind_2/"), pattern = "\\.log$")
files_structure_2 <- paste0(getwd(),"/clumpak_oryx_all_ind_2/",files_structure)
n_last_file_name <- 4   
n_first_line <- 23
df_likelihood <- as.data.frame(matrix(nrow = 8, ncol =20 ))
for(i in 1:length(files_structure_2)){
  file_name <- files_structure_2[i]
  k_run <-  as.numeric(substr(file_name, nchar(file_name) - n_last_file_name , nchar(file_name)- n_last_file_name))
  k_replicate <- as.numeric(str_match(file_name, "infile.run_\\s*(.*?)\\s*\\.")[2])
  likelihood <- as.character(unname(unlist(read.pattern(file =files_structure_2[i], pattern = "^Marginal Likelihood = .*"))))
  likelihood_2 <- as.numeric(substr(likelihood, n_first_line , nchar(likelihood)))
  df_likelihood[k_run,k_replicate] <- likelihood_2
}
df_likelihood_res <- rowMeans(df_likelihood)

ggplot()+
  geom_line(aes(x=1:8,y=df_likelihood_res),size=1) +
  geom_point(aes(x=1:8,y=df_likelihood_res),size=2,color="blue")+
  theme_bw(base_size = 14) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom")+
  xlab("K") +
  ylab("Marginal Likelihood")+
  scale_x_continuous(breaks = round(seq(1, 8, by = 1),1))
##################################
############# PCOA ###############
##################################
##################################
############## ALL ###############
##################################

colors_oryx <- c("darkgoldenrod1","darkslategray3","deeppink","goldenrod4","dodgerblue")
col2hex(colors_oryx)

ind_names_all <- oryx_gl$ind.names
ind_to_plot_all <- c("505","523","477","486","564","573","410","524","450","495","471","533","493","562","576","454","495","517","470","500","431","452","484","490")
loc_to_plot_all <- which(sapply(ind_names_all, function(x) x %in% ind_to_plot_all))
oryx_labels_PCOA_all <- oryx_gl
pcoa_oryx_2_all <- gl.pcoa(oryx_labels_PCOA_all)
oryx_labels_PCOA_all$ind.names <- ind_names_all
oryx_labels_PCOA_all$ind.names <- rep("" ,length(oryx_labels_PCOA_all$ind.names)) 
oryx_labels_PCOA_all$ind.names[as.numeric(loc_to_plot_all)] <- names(loc_to_plot_all)

gl.pcoa.plot(pcoa_oryx_2_all, oryx_labels_PCOA_all,labels="ind")

PCOA_test_2_all <- pcoa_oryx_2_all$scores
PCOA_test_2_all <-  as.data.frame(PCOA_test_2_all)
PCOA_test_2_all$pop <-  oryx_gl$pop

plot_ly(PCOA_test_2_all,x=~PC1,y=~PC2,z=~PC3,
        marker = list(size = 8),colors = colors_oryx)%>% 
  add_markers(color=~pop)
##################################
###### WITHOUT HISTORICAL ########
##################################
colors_wo_hist <- colors_oryx[3:5]

ind_nameswo_hist <- without_historical$ind.names
ind_to_plotwo_hist <- c("505","523","477","486","564","573","410","524","450","495","471","533","493","562","576","454","495","517","470","500","431","452","484","490")
loc_to_plotwo_hist <- which(sapply(ind_nameswo_hist, function(x) x %in% ind_to_plotwo_hist))
oryx_labels_PCOAwo_hist <- without_historical
pcoa_oryx_2wo_hist <- gl.pcoa(oryx_labels_PCOAwo_hist)
oryx_labels_PCOAwo_hist$ind.names <- ind_nameswo_hist
oryx_labels_PCOAwo_hist$ind.names <- rep("" ,length(oryx_labels_PCOAwo_hist$ind.names)) 
oryx_labels_PCOAwo_hist$ind.names[as.numeric(loc_to_plotwo_hist)] <- names(loc_to_plotwo_hist)

PCOA_test_2wo_hist <- pcoa_oryx_2wo_hist$scores
PCOA_test_2wo_hist <-  as.data.frame(PCOA_test_2wo_hist)
PCOA_test_2wo_hist$pop <-  without_historical$pop
plot_ly(PCOA_test_2wo_hist,x=~PC1,y=~PC2,z=~PC3,
        marker = list(size = 8),colors = colors_wo_hist)%>% 
  add_markers(color=~pop)
##################################
############ Oman_UAE ############
##################################
colors_oman_uae <- colors_oryx[4:5]
ind_namesoman_uae <- Oman_UAE$ind.names
ind_to_plotoman_uae <- c("505","523","477","486","564","573","410","524","450","495","471","533","493","562","576","454","495","517","470","500","431","452","484","490")
loc_to_plotoman_uae <- which(sapply(ind_namesoman_uae, function(x) x %in% ind_to_plotoman_uae))
oryx_labels_PCOAoman_uae <- Oman_UAE
pcoa_oryx_2oman_uae <- gl.pcoa(oryx_labels_PCOAoman_uae)
oryx_labels_PCOAoman_uae$ind.names <- ind_namesoman_uae
oryx_labels_PCOAoman_uae$ind.names <- rep("" ,length(oryx_labels_PCOAoman_uae$ind.names)) 
oryx_labels_PCOAoman_uae$ind.names[as.numeric(loc_to_plotoman_uae)] <- names(loc_to_plotoman_uae)

gl.pcoa.plot(pcoa_oryx_2oman_uae, oryx_labels_PCOAoman_uae)

PCOA_test_2oman_uae <- pcoa_oryx_2oman_uae$scores
PCOA_test_2oman_uae <-  as.data.frame(PCOA_test_2oman_uae)
PCOA_test_2oman_uae$pop <-  Oman_UAE$pop

plot_ly(PCOA_test_2oman_uae,x=~PC1,y=~PC2,z=~PC3,
        marker = list(size = 8),colors = colors_oman_uae)%>% 
  add_markers(color=~pop)
##################################
####### DIVERSITY ################
##################################
oryx_genind <- gl2gi(oryx_gl)
oryx_hierfstat <- genind2hierfstat(oryx_genind)
# q profile
oryx_diversity <- gl.report.diversity(oryx_gl)
# private alleles
oryx_private_alleles <- gl.report.pa(oryx_gl)
# FST
oryx_fst <- gl.fst.pop(oryx_gl)
#observed heterozygosity
Hobs <- gl.report.heterozygosity(oryx_gl)
# basic summary table for population genetic analyses
poppr_res <- poppr(oryx_genind,sample = 999)
hierfstat_res <- basic.stats(oryx_hierfstat)
observed_he <- colMeans(hierfstat_res$Ho,na.rm = T)
fis <- colMeans(hierfstat_res$Fis,na.rm = T)
expected_he <- colMeans(hierfstat_res$Hs,na.rm = T)

oryx_loci <- genind2loci(oryx_genind)
oryx_allelic_richness <- allelicrichness(oryx_loci)
oryx_rarefactionplot <- rarefactionplot(oryx_loci)
oryx_allelic.richness<- allelic.richness(oryx_hierfstat)
all_richness <- colMeans(oryx_allelic.richness$Ar)

##################################
############ Ne ##################
##################################
oryx_Ne <- write_genlight(oryx_vcf_radiator)
oryx_Ne <- gl.compliance.check(oryx_Ne)
oryx_Ne <- gl.recalc.metrics(oryx_Ne)
ploidy(oryx_Ne) <- 2
oryx_Ne <- gl.drop.ind(oryx_Ne, ind.list=individuals_to_remove)
# individual 454 was mislabelled and was changed from WWR_Oman to WWR-Mix
ind_454 <- which(oryx_Ne$ind.names=="OXL454")
oryx_Ne$pop[ind_454] <- as.factor("WWR_Mix")
oryx_Ne <- gl2gi(oryx_Ne)
all_together <- genomic_converter(oryx_genind,output = "hierfstat",path.folder = getwd())
##################################
############ HeatMAPS ############
##################################
##################################
###### WITHOUT HISTORICAL ########
##################################
without_historical$ind.names <- paste0(without_historical$ind.names,substring(without_historical$other$info$Gender,1,1))
colors_pops_all <- as.data.frame(without_historical$pop)
colors_pops_all$colors <- NA
colnames(colors_pops_all) <- c("pops","colors")
colors_pops_all[which(colors_pops_all$pops=="WWR Mix"),"colors"] <- col2hex(colors_oryx)[3]
colors_pops_all[which(colors_pops_all$pops=="WWR Oman"),"colors"] <- col2hex(colors_oryx)[4]
colors_pops_all[which(colors_pops_all$pops=="WWR UAE"),"colors"] <- col2hex(colors_oryx)[5]
oryx_related_matrix_all <- gl.grm(without_historical)
##################################
########### Oman_UAE #############
##################################
Oman_UAE$ind.names <- paste0(Oman_UAE$ind.names,substring(Oman_UAE$other$info$Gender,1,1))
colors_pops_uae_oman <- as.data.frame(Oman_UAE$pop)
colors_pops_uae_oman$colors <- NA
colnames(colors_pops_uae_oman) <- c("pops","colors")
colors_pops_uae_oman[which(colors_pops_uae_oman$pops=="WWR Oman"),"colors"] <- col2hex(colors_oryx)[4]
colors_pops_uae_oman[which(colors_pops_uae_oman$pops=="WWR UAE"),"colors"] <- col2hex(colors_oryx)[5]
oryx_related_matrix_uae_oman <- gl.grm(Oman_UAE)
##################################
############## MIX ###############
##################################
Mix$ind.names <- paste0(Mix$ind.names,substring(Mix$other$info$Gender,1,1))
colors_pops_mix <- as.data.frame(Mix$pop)
colors_pops_mix$colors <- NA
colnames(colors_pops_mix) <- c("pops","colors")
colors_pops_mix[which(colors_pops_mix$pops=="WWR Mix"),"colors"] <- col2hex(colors_oryx)[3]
oryx_related_matrix_mix <- gl.grm(Mix)
##################################
######## RELATEDNESS #############
##################################
oryx_related_matrix_all <- gl.grm(oryx_gl)
gl.grm.network(oryx_related_matrix_all,oryx_gl,node.label=T, method = "kk",node.label.size=1,node.size=10, node.label.color = "black",alpha = 0.006)
##################################
############# AMOVA ##############
##################################
##################################
############## ALL ###############
##################################
oryx_gl_2_all <- oryx_gl
oryx_genind_2_all <- gl2gi(oryx_gl_2_all)
# assign pop information to strata
strata(oryx_genind_2_all) <- data.frame(pop = pop(oryx_genind_2_all))
amova_all <- poppr.amova(oryx_genind_2_all,~pop,nperm = 10000)
##################################
###### WITHOUT HISTORICAL ########
##################################
oryx_gl_2_wo_hist <- without_historical
oryx_genind_2_wo_hist <- gl2gi(oryx_gl_2_wo_hist)
# assign pop information to strata
strata(oryx_genind_2_wo_hist) <- data.frame(pop = pop(oryx_genind_2_wo_hist))
amova_wo_hist <- poppr.amova(oryx_genind_2_wo_hist,~pop,nperm = 10000)
##################################
############ Oman_UAE ############
##################################
oryx_gl_2_uae_oman <- Oman_UAE
oryx_genind_2_uae_oman <- gl2gi(oryx_gl_2_uae_oman)
# assign pop information to strata
strata(oryx_genind_2_uae_oman) <- data.frame(pop = pop(oryx_genind_2_uae_oman))
amova_uae_oman <- poppr.amova(oryx_genind_2_uae_oman,~pop,nperm = 10000)

##################################
######### RELATEDNESS ############
##################################
##################################
######### SEQUOIA ################
##################################

LifeHistData_2 <- as.data.frame(oryx_gl$ind.names)
LifeHistData_2$Sex <- as.character(oryx_gl$other$info$Gender)
LifeHistData_2[LifeHistData_2=="Female"] <- 1
LifeHistData_2[LifeHistData_2=="Male"] <- 2
LifeHistData_2$Sex <- as.numeric(LifeHistData_2$Sex)

LifeHistData_2$BirthYear <- as.character(oryx_gl$other$info$Age)
LifeHistData_2[LifeHistData_2=="Old"] <- 2000
LifeHistData_2[LifeHistData_2=="Adult"] <- 2005
LifeHistData_2[LifeHistData_2=="Juvenile"] <- 2010
LifeHistData_2[LifeHistData_2=="Calf"] <- 2015
LifeHistData_2$BirthYear <- as.numeric(LifeHistData_2$BirthYear)

colnames(LifeHistData_2)<- c("ID","Sex","BirthYear")

oryx_sequoia <- as.matrix(as.data.frame(oryx_gl))
oryx_sequoia[is.na(oryx_sequoia)] <- "-9"
oryx_sequoia_2 <- cbind(oryx_gl$ind.names,oryx_sequoia)
oryx_sequoia_2 <- mapply(oryx_sequoia_2, FUN=as.numeric)
oryx_sequoia_2 <- matrix(data=oryx_sequoia_2, ncol=ncol(oryx_sequoia)+1, nrow=nrow(oryx_sequoia))
oryx_sequoia_2 <- as.data.frame(oryx_sequoia_2)
GenoM <- GenoConvert(InFile=oryx_sequoia_2, InFormat="seq",Missing="-9",IDcol=1)

GenoM_2 <- GenoM

sequoia_res <- sequoia(GenoM = GenoM_2,Plot = F,UseAge = F,MaxSibshipSize=20)
sequoia_res_2 <- GetMaybeRel(GenoM = GenoM_2,  LifeHistData =LifeHistData_2)

pedigree_res <- sequoia_res$Pedigree
pedigree_res$sex <- c(as.character(oryx_gl$other$info$Gender),rep("Female",3),rep("Male",3))
family <- pedigree(id=pedigree_res$id, dadid=pedigree_res$sire, momid=pedigree_res$dam,sex=pedigree_res$sex)
plot(family)

##################################
######### INBREEDING #############
##################################

inbreeding_oryx <- read.csv("inbreeding_oryx.csv")
Mix_inbreed <- inbreeding_oryx[which(inbreeding_oryx$Group=="WWR Mix"),]
UAE_inbreed <- inbreeding_oryx[which(inbreeding_oryx$Group=="WWR UAE"),]
Oman_inbreed <- inbreeding_oryx[which(inbreeding_oryx$Group=="WWR Oman"),]

res_inbreed <- as.data.frame(cbind(summary(Mix_inbreed)[,5:6],summary(UAE_inbreed)[,5:6],summary(Oman_inbreed)[,5:6]))
