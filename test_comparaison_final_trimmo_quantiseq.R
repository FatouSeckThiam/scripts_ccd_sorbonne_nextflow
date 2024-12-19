library(xCell)
library(readr)
library(edgeR)
library(immunedeconv)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(devtools)
library(DESeq2)
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)

#remarque = Immunocompétent -> IC
#Immuno Status -> Immune Status
#faire apparaitre la valeur de p pour HIV/IC neutrophils et faire descendre les étoiles (Karim) 
# HIV/IC POUMON VOIR LES MATRICES DA,S LES 2 FIGURES DE BAPTISTE(b CELLES ET NEUTROPHILS) A REVOIR !!!!!!!!
# Lib size enlever 200-SS, 063, 198
# ajouter une légende pour différencier les localisation poumons des autres dans le plot libsize 3p ( sans couleurs) et enlever l'ID en plus (n=20)

### papier =  > 
rm(list=ls())

##### TEST COMPARAISON FINAL Baptiste Quantiseq TRIMMOMATIC ######
p.workDir <- "/home/fatou/Documents/Baptiste_Poumon/Final_files_poumon/deconv_res_scripts_matrices/Trimmomatic_Quantiseq/"
#res_quantiseq_3p_trimmo_clinical_elispot = read.csv(paste0(p.workDir, "quantiseq_3p_poumon_trimmomatic_cell_fractions_infos_clinical.csv"),
#                                                  sep="\t", header = TRUE, as.is=TRUE)

### Enlever la colones Others
#res_quantiseq_3p_trimmo_clinical_elispot <- select(res_quantiseq_3p_trimmo_clinical_elispot, -Other)

#normalisation des données de sorte que la somme des types cellulaires pour chaque sample soit égale à 1 par patient
#res_quantiseq_3p_trimmo_clinical_elispot[, -1] <- res_quantiseq_3p_trimmo_clinical_elispot[, -1] / rowSums(res_quantiseq_3p_trimmo_clinical_elispot[, -1], na.rm = TRUE)

#write.table(res_quantiseq_3p_trimmo_clinical_elispot,"/home/fatou/Documents/Baptiste_Poumon/Final_files_poumon/deconv_res_scripts_matrices/Trimmomatic_Quantiseq/quantiseq_3p_poumon_trimmomatic_cell_fractions_normalised_infos_clinical.csv", row.names = F, quote = F, sep="\t")

#  matrice directement normalisée avec les nadir_cd4, res_ELISPOT, statut_immuno_nadir,	statut_immuno,localization
res_quantiseq_3p_trimmo_clinical_elispot = read.csv(paste0(p.workDir, "quantiseq_3p_poumon_trimmomatic_cell_fractions_normalised_infos_clinical.csv"),
                                                    sep="\t", header = TRUE, as.is=TRUE)

#Pour vérifier que toutes les lignes sont à 1
#res_quantiseq_3p_trimmo_clinical_elispot <- res_quantiseq_3p_trimmo_clinical_elispot %>% select(-c(nadir_cd4, res_ELISPOT, statut_immuno_nadir,	statut_immuno,localization))

#res_quantiseq_3p_trimmo_clinical_elispot = read.csv("/home/fatou/Téléchargements/quantiseq_3p_poumon_trimmomatic_cell_fractions_normalised_infos_clinical.csv",
                                                    #sep="\t", header = TRUE, as.is=TRUE)

# Transformer la matrice en long tout en conservant les colonnes sample, nadir_cd4, statut_immuno, 
res_quantiseq_3p_trimmo_clinical_elispot <- res_quantiseq_3p_trimmo_clinical_elispot %>%
  gather(cell_type, fraction, -Sample, -nadir_cd4, -res_ELISPOT, -statut_immuno_nadir, -statut_immuno, -localization)

dim(res_quantiseq_3p_trimmo_clinical_elispot)
### Créer des subsets
subset_vih_trimmo_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.01.193.ZJ", "FR.04.041.DJP", "FR.04.060.SJB","FR.04.062.GC",
                                                                                           "FR.04.063.FC",  "FR.04.127.FG", "FR.04.144.OH",  "FR.04.145.YK",  "FR.04.133.GB", "FR.01.090.SP", "FR.01.092.PA"))
#dim(subset_vih_trimmo_all)
subset_vih_trimmo_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.01.193.ZJ", "FR.04.060.SJB","FR.04.062.GC",
                                                                                              "FR.04.063.FC",  "FR.04.127.FG", "FR.04.144.OH",  "FR.04.145.YK",  "FR.04.133.GB"))

#dim(subset_vih_trimmo_poumon)
subset_ic_trimmo_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.01.162.PG", "FR.01.173.RP", "FR.01.174.PDRJ", "FR.01.177.BA", "FR.01.197.SC",
                                                                                          "FR.01.198.BL", "FR.04.196.MC", "FR.01.171.AMM", "FR.01.175.JF"))
#dim(subset_ic_trimmo_all)

subset_ic_trimmo_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.01.173.RP", "FR.01.174.PDRJ", "FR.01.177.BA", "FR.01.197.SC",
                                                                                             "FR.01.198.BL", "FR.04.196.MC", "FR.01.171.AMM"))
#dim(subset_ic_trimmo_poumon)

yes_elispot_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.01.173.RP", "FR.01.174.PDRJ", "FR.01.175.JF", "FR.01.162.PG", "FR.04.144.OH",
                                                                                     "FR.04.060.SJB", "FR.04.062.GC", "FR.04.063.FC" ))

no_elispot_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.04.196.MC", "FR.01.198.BL","FR.01.197.SC", "FR.01.177.BA", "FR.01.171.AMM",
                                                                                    "FR.04.133.GB", "FR.04.127.FG", "FR.04.041.DJP", "FR.01.193.ZJ", "FR.01.092.PA"))

yes_elispot_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample%in%c("FR.01.173.RP", "FR.01.174.PDRJ", "FR.04.144.OH",
                                                                                      "FR.04.060.SJB", "FR.04.062.GC", "FR.04.063.FC"))                                                                          

no_elispot_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.04.196.MC", "FR.01.198.BL","FR.01.197.SC", "FR.01.177.BA", "FR.01.171.AMM",
                                                                                       "FR.04.133.GB", "FR.04.127.FG", "FR.01.193.ZJ"))                                                                                 

HIV_nadir_haut_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in%c("FR.04.060.SJB","FR.04.062.GC", "FR.04.063.FC", "FR.01.090.SP", "FR.01.193.ZJ")) # FR-05-192-TD à vérifier car absent de la liste poumon

HIV_nadir_haut_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in%c("FR.04.060.SJB","FR.04.062.GC", "FR.04.063.FC", "FR.01.193.ZJ")) # FR-05-192-TD à vérifier car absent de la liste poumon

HIV_nadir_bas_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.04.145.YK", "FR.04.133.GB", "FR.04.127.FG", "FR.04.041.DJP", "FR.01.090.SP"))

HIV_nadir_bas_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in% c("FR.04.145.YK", "FR.04.133.GB", "FR.04.127.FG"))

autres_no_HIV_nadir_bas_poumon <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(Sample %in%c("FR.04.060.SJB","FR.04.062.GC", "FR.04.063.FC", "FR.01.193.ZJ",
                                                                                                   "FR.01.173.RP", "FR.01.174.PDRJ", "FR.01.177.BA",
                                                                                                   "FR.01.197.SC", "FR.01.198.BL", "FR.04.196.MC", "FR.01.171.AMM"))

autres_no_HIV_nadir_bas_all <- res_quantiseq_3p_trimmo_clinical_elispot %>% filter(!Sample %in% c("FR.04.145.YK", "FR.04.133.GB", "FR.04.127.FG", "FR.04.041.DJP", "FR.01.092.PA", "FR.04.144.OH"))                                                                                                                                                                  

### Test de comparaison wilcox
# ) L'ensemble des 20 patients: 
# 1) Répondeurs (175 / 174 / 173 / 162 / 144 / 063 / 062 / 060) vs non répondeurs (198 / 196 / 197 / 177 / 171 / 133 / 127 / 041 / 193 / 092)
elispot_all <- bind_rows(no_elispot_all, yes_elispot_all)  

#result <- compare_means(fraction ~ res_ELISPOT , group.by = "cell_type",  data = elispot_all,  method = "wilcox.test", p.adjust.method = "fdr")

# s'assurer que res_ELISPOT est bien de type factor 
if (!is.factor(elispot_all$res_ELISPOT)) {
  elispot_all$res_ELISPOT <- as.factor(elispot_all$res_ELISPOT)
}

# Filtrer les NA dans res_ELISPOT et obtenir les lignes uniques par patient
total_responses <- elispot_all %>%
  filter(!is.na(res_ELISPOT)) %>%  # Filtrer les NA
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(res_ELISPOT) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$res_ELISPOT, " n = ", total_responses$total ),
  total_responses$res_ELISPOT
)

p <- ggplot(elispot_all, aes(x = cell_type, y = fraction, fill = res_ELISPOT)) +
  geom_boxplot(outlier.shape = NA) +  # Enlever les points outliers
  theme_classic2() +
  labs(x = "Cell Type", y = "Cell fraction", fill = "ELISPOT-IFNγ (All samples)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("deeppink3", "yellowgreen"), label = label_map)

# Ajouter les résultats du test de Wilcoxon et #Ajouter les nombres de yes_response et no_response
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", 
                       label.y = max(elispot_all$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method = "fdr", 
                     label.y = max(elispot_all$fraction) * 1.2) 



#  Pour les patients issus uniquement d'un prélèvement poumon: 
# 2) Répondeurs (174 / 173 / 144 / 063 / 062 / 060) vs non répondeurs (198 / 196 / 197 / 177 / 171 / 133 / 127 / 193)

elispot_poumon <- bind_rows(no_elispot_poumon, yes_elispot_poumon)

total_responses <- elispot_poumon %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(res_ELISPOT) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$res_ELISPOT, " n = ", total_responses$total ),
  total_responses$res_ELISPOT
)

p <- ggplot(elispot_poumon, aes(x = cell_type, y = fraction, fill = res_ELISPOT)) +
  geom_boxplot(outlier.shape = NA) +  # Enlever les points outliers
  theme_classic2() +
  labs(x = "Cell Type", y = "Cell fraction", fill = "ELISPOT-IFNγ (Lung)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("deeppink3", "yellowgreen"), label = label_map)

# Ajouter les résultats du test de Wilcoxon et #Ajouter les nombres de yes_response et no_response
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", 
                       label.y = max(elispot_poumon$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method = "fdr", 
                     label.y = max(elispot_poumon$fraction) * 1.2) 

# L'ensemble des 20 patients: 
# 3) IC (198 / 196 / 197 / 177 / 175 / 174 / 173 / 171 / 162)
# vs HIV nadir haut (063 / 062 / 060 / 192 / 090) vs HIV nadir bas (145 / 133 / 127 / 041 / 092)
nadir_all <- bind_rows(HIV_nadir_haut_all, HIV_nadir_bas_all)
IC_nadir_all <- bind_rows(subset_ic_trimmo_all, nadir_all)

total_responses <- IC_nadir_all %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(statut_immuno_nadir) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$statut_immuno_nadir, " n = ", total_responses$total ),
  total_responses$statut_immuno_nadir
)


p <- ggplot(IC_nadir_all, aes(x = cell_type, y = fraction, fill = statut_immuno_nadir)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell type", y = "cell fraction", fill = "Immuno-Status-Nadir-CD4 (All samples)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("goldenrod2", "plum2", "turquoise2"),label = label_map)
  
  
# Ajouter les résultats du test de Wilcoxon et #Ajouter les nombres de yes_response et no_response
p + stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods = "fdr", 
                                             label.y = max(IC_nadir_all$fraction) * 1.1) 

# Calculer les p-values et les stocker dans un tableau
stat_test <- compare_means(fraction ~ statut_immuno_nadir, data = IC_nadir_all, 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")

# Afficher le tableau dans la console
print(stat_test)

                                                                                                                                                                              
# le test de Wilcox compare 2 à 2 les variable et vu qu'n a 3 variables, on a donc 2 comparaisons et donc 2 p-values par groupe 
#Et ces derniers se superposent dans le plot donc mieux afficher que le p.signif et pas les pvalues pour ne pas trop charger la figure
# p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", 
#                        label.y = max(IC_nadir_all$fraction) * 1.1) + 
#   stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method = "fdr", 
#                      label.y = max(IC_nadir_all$fraction) * 1.2) +  # Ajuster label.y pour p.signif au-dessus de p
#   # Ajuster label.y si nécessaire) +
#   annotate("text", x = Inf, y = Inf, label = paste("nadir_haut n =", total_responses$total[total_responses$statut_immuno_nadir == "HIV_with_high_CD4_nadir"],
#                                                    "nadir_bas n =", total_responses$total[total_responses$statut_immuno_nadir == "HIV_with_low_CD4_nadir"]),
#           

# 4) HIV nadir bas (145 / 133 / 127 / 041 / 092) 
# vs autres (198 / 196 / 197 / 177 / 175 / 174 / 173 / 171 / 162 / 063 / 062 / 060 / 192 / 090) # à rajouter dire à Baptiste 193 

nadir_bas_autres_all <- bind_rows(HIV_nadir_bas_all, autres_no_HIV_nadir_bas_all)

#la colonne nouveau_statut est créée, où les niveaux HIV-nadir_haut et IC sont fusionnés en un seul niveau appelé HIV-nadir_haut_IC.
nadir_bas_autres_all <- nadir_bas_autres_all %>%
  mutate(nouveau_statut = ifelse(statut_immuno_nadir %in% c("HIV_with_high_CD4_nadir", "Immunocompetent"),
                                 "Autres", statut_immuno_nadir))


total_responses <-nadir_bas_autres_all  %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(nouveau_statut) %>%
  summarise(total = n())


label_map <- setNames(
  paste0(total_responses$nouveau_statut, " (n = ", total_responses$total ),
  total_responses$nouveau_statut
)

#Enlever les valeurs nulles
#nadir_bas_autres_all <- nadir_bas_autres_all %>% filter(nouveau_statut != "")

p <- ggplot(nadir_bas_autres_all, aes(x = cell_type, y = fraction, fill = nouveau_statut )) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell type", y = "Cell fraction", fill = "Immuno-Status-Nadir-CD4 (All samples)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("royalblue", "plum2"), label = label_map)
  
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", # ne pas enlever dans la figure VIH/IC poumons
                       label.y = max(nadir_bas_autres_all$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method = "fdr", 
                     label.y = max(nadir_bas_autres_all$fraction) * 1.2) 

stat_test <- compare_means(fraction ~ nouveau_statut, data = nadir_bas_autres_all, 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")

# Afficher le tableau dans la console
print(stat_test)

### Pour les patients issus uniquement d'un prélèvement poumon: 

### 5) - IC (198 / 196 / 197 / 177 / 174 / 173 / 171) vs HIV nadir haut (063 / 062 / 060 / 192) vs HIV nadir bas (145 / 133 / 127)
nadir_poumon <- bind_rows(HIV_nadir_haut_poumon, HIV_nadir_bas_poumon)
IC_nadir_poumon <- bind_rows(subset_ic_trimmo_poumon, nadir_poumon)

# Important juste pour vérifier qu'on compare les bons groupes et la manière de comparaison.
#Pour cette comparaison de 3, le test wilcox compare 2 à 2, donc 3 combinaisons par cellule. 
# Si pour une cellule, les 3 combinaisons sont ns, donc cest ns dans la figure
# et si une combinaison donne un résultat significatif parmis les 3, donc le plot sort * et ns en mm temps. 
# et c'est le cas de ce cas de figure avec les Neutrophiles: 
# 13 Neutrophils    fraction Immunocompetent         HIV_with_high_CD4_nadir 0.106      1 0.106    ns       Wilcoxon
# 14 Neutrophils    fraction Immunocompetent         HIV_with_low_CD4_nadir  0.0456     1 0.046    *        Wilcoxon
# 15 Neutrophils    fraction HIV_with_high_CD4_nadir HIV_with_low_CD4_nadir  0.505      1 0.505    ns       Wilcoxon
#result <- compare_means(fraction ~ statut_immuno_nadir  , group.by = "cell_type",  data = IC_nadir_poumon,  method = "wilcox.test", p.adjust.method = "fdr")

total_responses <- IC_nadir_poumon %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(statut_immuno_nadir) %>%
  summarise(total = n())

label_map <- setNames(
  paste0(total_responses$statut_immuno_nadir, " n = ", total_responses$total),
  total_responses$statut_immuno_nadir
)


p <- ggplot(IC_nadir_poumon, aes(x = cell_type, y = fraction, fill = statut_immuno_nadir)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell type", y = "cell fraction", fill = "Immune Status Nadir CD4 (Only Lung)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("goldenrod2", "plum2", "turquoise2"), label = label_map)


# Ajouter les résultats du test de Wilcoxon et #Ajouter les nombres de yes_response et no_response
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", # ne pas enlever dans la figure VIH/IC poumons
                       label.y = max(nadir_bas_autres_all$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method = "fdr", 
                     label.y = max(nadir_bas_autres_all$fraction) * 1.2)


# Calculer les p-values et les stocker dans un tableau
stat_test <- compare_means(fraction ~ statut_immuno_nadir, data = IC_nadir_poumon , 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")


# Sauvegarder les résultats dans un fichier CSV
write.csv(stat_test, "wilcoxon_test_results_nadir_immuno_poumon.csv", row.names = FALSE)


### 6) - HIV nadir bas (145 / 133 / 127) vs autres (198 / 196 / 197 / 177 / 174 / 173 / 171 / 063 / 062 / 060 / 192)
nadir_bas_autres_poumon <- bind_rows(HIV_nadir_bas_poumon, autres_no_HIV_nadir_bas_poumon)

#la colonne nouveau_statut est créée, où les niveaux HIV-nadir_haut et IC sont fusionnés en un seul niveau appelé HIV-nadir_haut_IC.
nadir_bas_autres_poumon <- nadir_bas_autres_poumon %>%
  mutate(nouveau_statut = ifelse(statut_immuno_nadir %in% c("HIV_with_high_CD4_nadir", "Immunocompetent"),
                                 "Autres", statut_immuno_nadir))


total_responses <-nadir_bas_autres_poumon  %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(nouveau_statut) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$nouveau_statut, " n = ", total_responses$total),
  total_responses$nouveau_statut
)

#Enlever les valeurs nulles
#nadir_bas_autres_all <- nadir_bas_autres_all %>% filter(nouveau_statut != "")

p <- ggplot(nadir_bas_autres_poumon, aes(x = cell_type, y = fraction, fill = nouveau_statut )) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell type", y = "Cell fraction", fill = "Immuno-Status-Nadir-CD4 (Only Lung)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("royalblue", "plum2"), label=label_map)

# Ajouter les résultats du test de Wilcoxon
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", 
                       label.y = max(nadir_bas_autres_poumon$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods = "fdr",
                     label.y = max(nadir_bas_autres_poumon$fraction) * 1.2)


# Calculer les p-values et les stocker dans un tableau
stat_test <- compare_means(fraction ~ nouveau_statut, data = nadir_bas_autres_poumon, 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")

### 7) comparaison HIV/IC all
vih_ic_all <- bind_rows(subset_vih_trimmo_all, subset_ic_trimmo_all)

# s'assurer que res_ELISPOT est bien de type factor 
if (!is.factor(vih_ic_all$statut_immuno)) {
  vih_ic_all$statut_immuno <- as.factor(vih_ic_all$statut_immuno)
}

# ÉLiminer l'espace qui s'est insérée dans la colonne statut_immuno variable HIV 
# Supprimer les espaces indésirables dans la colonne statut_immuno pour éliminer la répétition de HIV doublé 
vih_ic_all <- vih_ic_all %>%
  mutate(statut_immuno = str_trim(statut_immuno))

total_responses <- vih_ic_all %>%
  distinct(Sample, .keep_all = TRUE) %>%  # Obtenir les lignes uniques par patient
  group_by(statut_immuno) %>%
  summarise(total = n())


# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$statut_immuno, " n = ", total_responses$total),
  total_responses$statut_immuno
)


p <- ggplot(vih_ic_all, aes(x = cell_type, y = fraction, fill = statut_immuno)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell Type", y = "Cell fraction", fill = "VIH vs IC (All samples)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("lightsteelblue", "turquoise2"), label=label_map)

p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.methods = "fdr", 
                       label.y = max(vih_ic_all$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods = "fdr",
                     label.y = max(vih_ic_all$fraction) * 1.2)

# Calculer les p-values et les stocker dans un tableau
stat_test <- compare_means(fraction ~ statut_immuno, data =vih_ic_all, 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")

### 8-1) comparaison HIV/IC Poumon
vih_ic_poumon <- bind_rows(subset_vih_trimmo_poumon, subset_ic_trimmo_poumon)

# s'assurer que statut_immuno est bien de type factor 
if (!is.factor(vih_ic_poumon$statut_immuno)) {
  vih_ic_poumon$statut_immuno <- as.factor(vih_ic_poumon$statut_immuno)
}

# ÉLiminer l'espace qui s'est insérée dans la colonne statut_immuno variable HIV 
# Supprimer les espaces indésirables dans la colonne statut_immuno pour éliminer la répétition de HIV doublé 
vih_ic_poumon <- vih_ic_poumon %>%
  mutate(statut_immuno = str_trim(statut_immuno))

# Calculer le nombre total de réponses par statut immunologique
total_responses <- vih_ic_poumon %>%
  distinct(Sample, .keep_all = TRUE) %>%
  group_by(statut_immuno) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$statut_immuno, " n = ", total_responses$total),
  total_responses$statut_immuno
)

# Modifier le plot pour utiliser les étiquettes personnalisées dans la légende
p <- ggplot(vih_ic_poumon, aes(x = cell_type, y = fraction, fill = statut_immuno)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell Type", y = "Cell fraction", fill = "VIH vs IC (Only Lung)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("HIV" = "lightsteelblue", "Immunocompetent" = "turquoise2"),
                    labels = label_map)


### 8-2) comparaison HIV/IC Poumon matrice cybersortx de KARIM (/home/fatou/Téléchargement/CIBERSORTx_Job16_Results.txt )
rm(list=ls())
res_3p_cybersort_karim = read.csv("/home/fatou/Téléchargements/CIBERSORTx_Job16_Results.txt",
                                                    sep="\t", header = TRUE, as.is=TRUE)


# Sélectionner que les patients poumon IC
subset_ic_poumon_karim <- res_3p_cybersort_karim %>% filter(Mixture %in% c("FR01173RP", "FR01174PDRJ", "FR01177BA", "FR01197SC",
                                                                                             "FR04198BL", "FR04196MC", "FR01171AMM"))

# Sélectionner que les patients poumon VIH
subset_vih_poumon_karim <- res_3p_cybersort_karim %>% filter(Mixture %in% c("FR01193ZJ", "FR04060SJB","FR04062GC", "FR04063FC",  "FR04127FG", "FR04144OH",  "FR04145YK",  "FR04133GB"))
                                                                                              
# Regrouper les 2 matrices
subset_3p_cybersort_karim_vih_ic_poumon <- bind_rows(subset_vih_poumon_karim, subset_ic_poumon_karim)

# Sélectionner dans la matrice que les cellules immunitaires et les infos cliniques en enelevant les stats de cybersort x
subset_3p_cybersort_karim_vih_ic_poumon <- subset_3p_cybersort_karim_vih_ic_poumon %>%select(-c("P.value", "Correlation", "RMSE"))

# Transformer la matrice en long tout en conservant les colonnes sample, nadir_cd4, statut_immuno
subset_3p_cybersort_karim_vih_ic_poumon <- subset_3p_cybersort_karim_vih_ic_poumon %>%
  gather(cell_type, fraction, -Mixture, -nadir_cd4, -res_ELISPOT, -statut_immuno_nadir, -statut_immuno, -localization)


# s'assurer que statut_immuno est bien de type factor 
if (!is.factor(subset_3p_cybersort_karim_vih_ic_poumon$statut_immuno)) {
  subset_3p_cybersort_karim_vih_ic_poumon$statut_immuno <- as.factor(subset_3p_cybersort_karim_vih_ic_poumon$statut_immuno)
}

# ÉLiminer l'espace qui s'est insérée dans la colonne statut_immuno variable HIV 
# Supprimer les espaces indésirables dans la colonne statut_immuno pour éliminer la répétition de HIV doublé 
subset_3p_cybersort_karim_vih_ic_poumon <- subset_3p_cybersort_karim_vih_ic_poumon %>%
  mutate(statut_immuno = str_trim(statut_immuno))

# Calculer le nombre total de réponses par statut immunologique
total_responses <- subset_3p_cybersort_karim_vih_ic_poumon %>%
  distinct(Mixture, .keep_all = TRUE) %>%
  group_by(statut_immuno) %>%
  summarise(total = n())

# Créer une étiquette pour chaque statut immunologique avec le total correspondant
label_map <- setNames(
  paste0(total_responses$statut_immuno, " n = ", total_responses$total),
  total_responses$statut_immuno
)

# Modifier le plot pour utiliser les étiquettes personnalisées dans la légende
p <- ggplot(subset_3p_cybersort_karim_vih_ic_poumon, aes(x = cell_type, y = fraction, fill = statut_immuno)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic2() +
  labs(x = "Cell Type", y = "Cell fraction", fill = "VIH vs IC (Only Lung)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("HIV" = "lightsteelblue", "Immunocompetent" = "turquoise2"),
                    labels = label_map)


# Ajouter les tests de comparaison des moyennes
p + stat_compare_means(method = "wilcox.test", label = "p", p.adjust.method = "fdr", 
                       label.y = max(subset_3p_cybersort_karim_vih_ic_poumon$fraction) * 1.1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods = "fdr",
                     label.y = max(subset_3p_cybersort_karim_vih_ic_poumon$fraction) * 1.2)

# Calculer les p-values et les stocker dans un tableau
stat_test <- compare_means(fraction ~ statut_immuno, data =subset_3p_cybersort_karim_vih_ic_poumon, 
                           group.by = "cell_type", method = "wilcox.test", 
                           p.adjust.method = "fdr")

# Sauvegarder les résultats dans un fichier CSV
p.workDir <- "/home/fatou/Documents/Baptiste_Poumon/Final_files_poumon/deconv_res_scripts_matrices/Trimmomatic_Quantiseq/Plots_comparasion_Baptiste/plots_colors_Baptiste/"
write.csv(stat_test, paste0(p.workDir, "wilcoxon_test_results_IC_HIV_poumon.csv"), row.names = FALSE)

### Boxplot QC depth wes tumour
p.workDir <- "/home/fatou/Documents/"
mean_depth_wes_poumon = read.csv(paste0(p.workDir, "mean_depth_tumour_poumon.csv"),
                                                    sep="\t", header = TRUE, as.is=TRUE)

ggplot(data = mean_depth_wes_poumon, aes(x=WES, y=mean_depth_cov_Tumour_poumon, fill=WES)) +
  stat_boxplot(geom = "errorbar", # Boxplot with error bars 
               width = 0.2) +
  geom_boxplot(colour = "black", # Colors
               alpha = 0.9, outlier.colour = "red") +
  scale_y_continuous(name = "mean depth coverage") +  # Continuous variable label
  #scale_x_discrete(name = "Feed") +      # Group label
  ggtitle("mean depth coverage WES") + # Plot title
  theme(axis.line = element_line(colour = "black", # Theme customization
                                 size = 0.2)) + 
  scale_fill_manual(values = c("darkolivegreen", "brown3"))  





