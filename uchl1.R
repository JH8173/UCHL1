setwd('/home/jun/melanoma/')

#####################################################################################

data <- read.csv('Table_S1D.csv', na.strings = c("-", 'n/a', '[Not Available]', '[ERROR]',
                                                 '[Not Applicable]'))
summary(data)
data$Name <- as.character(data$Name)
data$UV.signature
data$MUTATIONSUBTYPES

data$Name[duplicated(gsub('-0[16]{1}$', '', data$Name))]

data <- data[data$Name != "TCGA-ER-A19T-06",]
data <- data[data$Name != "TCGA-ER-A2NF-06",]

ID_BRAF <- data$Name[data$MUTATIONSUBTYPES == 'BRAF_Hotspot_Mutants']
ID_RAS <- data$Name[data$MUTATIONSUBTYPES == 'RAS_Hotspot_Mutants']
ID_NF1 <- data$Name[data$MUTATIONSUBTYPES == 'NF1_Any_Mutants']
ID_TW_UV <- data$Name[data$MUTATIONSUBTYPES == 'Triple_WT' & 
                        data$UV.signature == 'UV signature']
ID_TW_NUV <- data$Name[data$MUTATIONSUBTYPES == 'Triple_WT' & 
                         data$UV.signature == 'not UV']
ID_UV <- data$Name[data$UV.signature == 'UV signature']
ID_NUV <- data$Name[data$UV.signature == 'not UV']

ID_BRAF <- ID_BRAF[!is.na(ID_BRAF)]
ID_RAS <- ID_RAS[!is.na(ID_RAS)]
ID_NF1 <- ID_NF1[!is.na(ID_NF1)]



###############################################################

library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
cs <- getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[78,1]

caselist <- getCaseLists(mycgds,mycancerstudy)

mycaselist = getCaseLists(mycgds,mycancerstudy)[7,1]

# Get available genetic profiles
geneticprofile = getGeneticProfiles(mycgds,mycancerstudy)

mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[2,1]

# Get data slices for a specified list of genes, genetic profile and case list
uchl1 <- getProfileData(mycgds,('UCHL1'),mygeneticprofile,mycaselist)


hist(uchl1$UCHL1, 
     breaks = 1000, 
     xlim = c(0,1000),
     plot = TRUE
     )

log(200, base = 10)
hist(log(uchl1$UCHL1))

##### http://www.gtexportal.org/home/gene/UCHL1 ##########

ID <- gsub('\\.', '-', (gsub('\\.0.', '', rownames(uchl1))))
df_uchl1 <- data.frame(ID = ID, UCHL1 = uchl1$UCHL1, 
                       log_UCHL1 = log(uchl1$UCHL1+1, base = 10),
                       UCHL1_G = factor(log(uchl1$UCHL1, base = 10) > 2.5)
                       )
colnames(data)
data$ID <- gsub('-0[0-9]{1}', '', data$Name)
data_uchl1 <- merge(df_uchl1, data, by = 'ID', all.x = TRUE)


colnames(data_uchl1)
boxplot(log_UCHL1 ~ UV.signature, data_uchl1)
t.test(log_UCHL1 ~ UV.signature, data_uchl1)

boxplot(log_UCHL1 ~ MUTATIONSUBTYPES, data_uchl1)
summary(aov(log_UCHL1 ~ MUTATIONSUBTYPES, data_uchl1))

boxplot(log_UCHL1 ~ ALL_PRIMARY_VS_METASTATIC, data_uchl1)
boxplot(log_UCHL1 ~ REGIONAL_VS_PRIMARY, data_uchl1)
boxplot(log_UCHL1 ~ CURATED_TCGA_SPECIMEN_SITE, data_uchl1)

a <- aov(log_UCHL1 ~ CURATED_TCGA_SPECIMEN_SITE, data_uchl1)
summary(a)
summary(data_uchl1)
