library(CancerInSilico)
library(methods)

## Full run examples

modDefault <- inSilicoCellModel(10,10,0.1)
modLongRun <- inSilicoCellModel(5,100,0.1)
modLargeRun <- inSilicoCellModel(1000,1,0.1)
modHighDensity <- inSilicoCellModel(100,10,0.3)

c1 <- new('CellType', name='DEFAULT', size=1, minCycle=48,
    cycleLength=function() 48)

c2 <- new('CellType', name='DOUBLE_SIZE', size=2, minCycle=48,
    cycleLength=function() 48)

d1 <- new('Drug', name='NO_EFFECT', timeAdded=0, cycleLengthEffect=
    function(type,len) len)

d2 <- new('Drug', name='HALF_CYCLE_LENGTH', timeAdded=0, cycleLengthEffect=
    function(type,len) len/2)

d3 <- new('Drug', name='HALF_DEFAULT_TYPE', timeAdded=0, cycleLengthEffect=
    function(type,len) ifelse(type=='DEFAULT', len/2, len))

d4 <- new('Drug', name='ADD_LATE', timeAdded=6, cycleLengthEffect=
    function(type,len) len)

modDrugs <- new('DrasdoHohmeModel', initialNum=100, runTime=4, density=0.4,
    randSeed=0, outputIncrement=4, recordIncrement=0.1, syncCycle=FALSE,
    cellTypes=c(c1,c2), cellTypeInitFreq=c(0.5,0.5), drugs=c(d1,d2,d3,d4))

save(modDefault, modLongRun, modLargeRun, modHighDensity, modCellTypes,
    modDrugs, file = 'SampleModels.RData')

## Pathway examples

# gene names

geneNamesGtoS <- c('AK2','ANP32E','ASF1A','ASF1B','ATAD2','AURKA','AURKB','BARD1','BIRC5','BRCA1','BRCA2','BRMS1L','BUB1B','CBX5','CCNB2','CCNE1','CDC20','CDC25A','CDC25B','CDCA3','CDCA8','CDK4','CDKN1A','CDKN1B','CDKN2A','CDKN2C','CDKN3','CENPE','CENPM','CHEK1','CHEK2','CIT','CKS1B','CKS2','CSE1L','CTCF','DCK','DCLRE1B','DCTPP1','DEK','DEPDC1','DIAPH3','DLGAP5','DNMT1','DONSON','DSCC1','DUT','E2F8','EED','EIF2S1','ESPL1','EXOSC8','EZH2','GINS1','GINS3','GINS4','GSPT1','H2AFX','H2AFZ','HELLS','HMGA1','HMGB2','HMGB3','HMMR','HN1','HNRNPD','HUS1','ILF3','ING3','IPO7','KIF18B','KIF22','KIF2C','KIF4A','KPNA2','LBR','LIG1','LMNB1','LYAR','MAD2L1','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','MELK','MKI67','MLH1','MRE11A','MSH2','MTHFD2','MXD3','MYBL2','MYC','NAP1L1','NASP','NBN','NCAPD2','NME1','NOLC1','NOP56','NUDT21','NUP107','NUP153','NUP205','PA2G4','PAICS','PAN2','PCNA','PDS5B','PHF5A','PLK1','PLK4','PMS2','PNN','POLA2','POLD1','POLD2','POLD3','POLE','POLE4','POP7','PPM1D','PPP1R8','PRDX4','PRIM2','PRKDC','PRPS1','PSIP1','PSMC3IP','PTTG1','RACGAP1','RAD1','RAD21','RAD50','RAD51AP1','RAD51C','RAN','RANBP1','RBBP7','RFC1','RFC2','RFC3','RNASEH2A','RPA1','RPA2','RPA3','RQCD1','RRM2','SHMT1','SLBP','SMC1A','SMC3','SMC4','SMC6','SNRPB','SPAG5','SPC24','SPC25','SSRP1','STAG1','STMN1','SUV39H1','SYNCRIP','TACC3','TBRG4','TCF19','TFRC','TIMELESS','TIPIN','TK1','TMPO','TOP2A','TP53','TRIP13','TUBB','TUBG1','UBE2S','UBE2T','UBR7','UNG','USP1','WDR90','WEE1','XPO1','XRCC6','ZW10','APAF1','ATM','CASP7','CCNA2','CCND3','CCNE2','CDC6','CDK2','CEBPA','CES1','CES2','CES3','CREBBP','DHFR','E2F1','E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','EP300','HBP1','HDAC1','HIC1','KAT2A','KAT2B','MCL1','PLAU','POLA1','PRMT5','RB1','RBBP4','RBBP8','RBL1','RBL2','RRM1','RYBP','SERPINE1','SIRT1','SMARCA2','SP1','SULT2A1','TFDP1','TFDP2','TFE3','TOPBP1','TP73','TRIM28','TRRAP','TYMS','UXT','WASF1','XRCC1','YY1','ABCC10','CCND1','CSF1','CSF2','FSHR','GGH','IGF1','MAT2A','PBK','PCSK6','SLC44A1','SP3','TERT','AIM1','BCL2L11','BCL2','BIK','CDKN2D','CYC1','DBF4','DIRAS3','DLK1','DNMT3A','DUSP4','FGFR1','FGFR2','FLG','GADD45B','HIST1H2AG','HIST1H2AK','HIST1H2AL','HIST1H2AM','HSPA5','IRF3','KCNH1','MAP3K5','MCPH1','MERTK','MYCN','NOX4','NRP1','NUSAP1','POLE2','POMC','PTPRO','SIVA1','TAPBP','TNFSF11','JAM2')
geneNamesGtoM <- c('ABL1','AKT1','BIRC2','BIRC3','CCNA2','CCND1','CCND2','CCND3','CCNE1','CDK2','CDK4','CDK6','CDKN1A','CDKN1B','CDKN2A','CREBBP','CSF2','E2F1','E2F2','E2F3','E2F4','EP300','FKBP1A','GRB2','GSK3B','HDAC1','JUN','MAP2K3','MAP3K7','MAPK14','MAPK9','MAPKAP1','MDM2','NFKB1','PML','PPP1CA','PPP1R15A','PPP2CA','PPP2CB','PRKDC','RAF1','RB1','RELA','RHOA','RIPK1','RPS6KB1','SFN','SKIL','SMAD2','SMAD3','SMAD4','SMAD7','SMURF1','SMURF2','SPTBN1','SQSTM1','SUV39H1','TFDP1','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','TNF','TNFAIP3','TNFRSF1A','TRAF1','TRAF2','XIAP','YWHAB','YWHAE','YWHAG','YWHAH','YWHAQ','YWHAZ','ZFYVE16','ZFYVE9')
geneNamesProx <- c('AMOT','AMOTL1','AMOTL2','CASP3','DVL2','LATS1','LATS2','NPHP4','SAV1','STK3','STK4','TJP1','TJP2','WWC1','WWTR1','YAP1','YWHAB','YWHAE','TNNT2','ACTA1','ACTA2','ADRA1A','ANKRD1','CKM','CSH2','CTGF','DES','FOXO3','GNG3','HGF','IFITM3','IL6','MSLN','MYF5','MYF6','MYH6','MYH7','MYL2','MYOCD','S100B','SLC2A4','TGFA','TNC','TNNC1','TPM1','CIITA','RCAN2','RCN2','PAX3','PCYT1A','CSH1','HSD3B1','APC2','YWHAQ','RASSF1','FZD10','FRMD6','WTIP','CSNK1D','CSNK1E','CTNNA1','CTNNA2','CTNNB1','GDF7','RASSF6','DLG1','AFP','DLG2','DLG3','DLG4','DVL1','DVL3','FGF1','FBXW11','CRB1','SCRIB','FZD2','TAZ','AMH','MIF','BBC3','GLI2','CRB2','CTNNA3','GSK3B','APC','BIRC2','BIRC5','ID1','ID2','BMP8A','ITGB2','AREG','AR','GDF6','LLGL2','LLGL1','SMAD1','SMAD2','SMAD3','SMAD4','SMAD7','MYC','NF2','SERPINE1','PARD6A','LEF1','WNT16','WNT4','PPP1CA','PPP1CB','PPP1CC','PPP2CA','PPP2CB','PPP2R1A','PPP2R1B','PPP2R2A','PPP2R2B','PPP2R2C','PRKCI','PPP2R2D','PRKCZ','PARD3','ASIP','CCND1','ACTB','MPP5','BMP2','BMP4','BMP5','BMP6','BMP7','BMP8B','BMPR1A','BMPR1B','BMPR2','SNAI2','SOX2','TCF7','TCF7L2','TCF4','TEAD1','TEAD4','TEAD3','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','ACTG1','TP53BP2','TP73','WNT1','WNT2','WNT3','WNT5A','WNT6','WNT7A','WNT7B','WNT8A','WNT8B','WNT10B','WNT11','WNT2B','WNT9A','WNT9B','YWHAG','YWHAH','YWHAZ','FZD5','FZD3','FRMD1','WNT10A','WNT5B','GDF5','AXIN1','AXIN2','FZD1','FZD4','FZD6','FZD7','FZD8','FZD9','TCF7L1','TCF3','PARD6G','PARD6B','TEAD2','NKD1','CCND2','BTRC','CCND3','WNT3A','LIMD1','CDH1')
geneNamesGrowth <- c('A2M','ACAT1','AGT','BACE1','C2','C7','CASP4','CCL19','CCL2','CCL5','CCND1','CD40','CD86','CDKN1A','CEBPD','CFB','CIITA','CSN2','CXCL10','CXCL9','FASN','FCER2','FCGR1A','FCGR1C','FCGRT','FOS','GLS','HAMP','HAS2','HBG1','ICAM1','IDO1','IFNG','IL12A','IL1B','IL21','IL2RA','IL6ST','IRF1','IRF8','ISG15','IVL','KLF4','LCN2','LPL','LY6E','LY96','MSR1','MUC1','MUC5B','MVP','MYC','MYD88','NOS2','NOX1','NPSR1','NR1H4','OPRM1','PEMT','PIM1','PPARGC1B','PSMB9','PTGS2','S100P','SCGB3A1','SOAT1','SOCS3','SPIN3','TAP1','TBX21','TFAP2A','TREX1','TRH','TYMP','VIP','WARS','XAF1','ZFP36','CEBPB','GBP1','IL10','IL18BP','PTN','WARS2','CDKN1B','SERPINA3','AHR','BIRC5','CCR5','CD274','CD46','CISH','CLGN','CRP','CTGF','CYP19A1','DMBT1','DUSP1','EGR1','EPO','FAAH','FGB','FGG','FGL2','FOXM1','GADD45B','HGF','HMGA2','HMOX1','HP','HSP90AA1','IL17A','IL1R1','IL2RG','JAK3','LBP','MCL1','MMP1','MMP2','MMP7','MPP2','NDN','NOS3','PFN1','PHB','PML','POMC','PRF1','REG1A','ROR1','SALL4','SERPINA1','SLC9A3','SPI1','SPIN1','SPIN2A','STAT3','TIMP1','TP53','TP63','TSPO','TWIST1','VEGFA','VEGFC','VIM','S1PR1','ACO1','BCL2L1','CCND2','FGF21','GNLY','IRF4','KIR3DL1','MUM1','ONECUT1','PAX5','SKP2','SLC10A1','SLC30A2','TLR2','ADH1C','C4A','CD2','CUX2','CYP4A11','IGF1','IGFALS','OSM','PIP','SOCS2','ACHE','CCT8','CD28','DDIT4','DPP3','EGR2','EIF2AK1','FUT1','FUT4','INSIG2','JAM2','KATNB1','KCNK5','LPA','MYLK','NR1H2','NR4A1','PDGFB','PRKACA','PRKCA','PSEN1','PYGO2','SOD1','SPAST','SRF','SULT1A1','SULT1A3','TAGLN','TLR4','TLR9','TNF','VIL1','ZNF225','ACTB','ELK1','ABCC1','ADIPOR1','ATF3','ATF4','BMI1','BOLA2','BRD7','CAD','CDK4','EGFR','EIF4E','FHIT','GATA4','ITGA7','LGALS1','MCM2','MCM7','MXD4','MYBBP1A','MYCT1','NHP2','NOP16','NR5A1','PDCD10','POLA2','PROM1','RPL10A','RPS19','RPS7','SERPINI1','SHBG','SKP1','SPP1','TERT','TFDP1','TIMP2','CDC25A','CXCR4','DDX18','HMGA1','LAP3','ODC1','PTMA','YBX1','BIRC7','MDM2','NES','PTDSS1','PTK2','SLC1A2','ABCB1','AHRR','ANG','APP','AR','ATP12A','ATP6AP2','BAK1','BCL2L10','BDKRB1','BDKRB2','BLK','BLNK','BMP2','BTRC','CCL20','CCL22','CD38','CD80','CEBPA','CES2','CHKA','CIDEA','COL1A2','CPB2','CREB1','CSF1','CSF2','CSNK1A1','CXCL1','CXCL2','CXCL3','CYP1A1','DIABLO','EBI3','EDN1','F3','F8','FASLG','FAS','FBP1','FGFR4','GCLC','GCLM','GH1','GIP','GNAI2','HBE1','HDAC8','HK3','HLA-C','HRAS','HSPA1A','IFNB1','IKBKB','IKZF3','IL12B','IL23A','IL2','IL32','IL6','ILK','ITM2B','KDR','KLK3','LIPG','LTA','MAT2A','MMP9','MYOZ2','NFKB1','NFKBIB','NOD2','NOLC1','NPY1R','NR3C1','NRF1','PDCD4','PIK3AP1','PLAUR','PLCD1','PTGES','S100A6','SAA1','SAA2','SCN5A','SELE','SELP','SIN3A','SIRT1','SLC2A4','SLC3A2','SOD2','SPR','ST8SIA1','TACR1','TANK','TF','TLR7','TNFAIP3','TNFSF13','TRAF2','UBE2D3','UBE2M','UGT1A1','UPK1B','VCAM1','ALDH1A1','COL11A1','IGF1R','OLFM4','ABCG2','ALOX12','APOE','BCL2A1','BDNF','CRMP1','HEPH','ING2','PTPN6','REL','SLC22A4','STAT4','TNFRSF10B','CD40LG','GUCY1A3','GUCY1B3','IFNGR2','MICA','NR4A2','AGER','BCL2','CREB3','CTCF','CTSB','ELF3','HIF1A','HLA-G','IL5','NQO1','OLR1','GRK5','IER3','PRG1','TNFSF15','ADAMTS5','ADORA2A','AGTR1','ALCAM','ALOX5AP','B2M','BCL10','BLVRA','BRD3','CDK2','CDX2','CHI3L1','COL4A3BP','CX3CL1','CXCL11','CYBA','DIO2','FLT3','FN1','FSTL3','FTH1','G6PC','GADD45A','GNRH2','HLA-A','IKBKE','IL18','IL1RN','IL27','IRAK3','LTB','LZTS2','MAN2B1','MGMT','MIR155HG','MUC2','MUC7','NAF1','NAMPT','NFKBIA','NGB','NLRP3','OPTN','OXTR','P2RY2','PASK','PCYT2','PDGFD','PHEX','PI3','PIGR','PLA2G2A','PLK1','POLK','POLQ','PPP1R13L','PRDX5','PRDX6','PTEN','PTGDS','RGS4','SCNN1A','SDC4','SERPINE1','SLC6A6','SNAI1','SOX9','SRA1','STK39','TAC1','TEP1','TFF1','TNFRSF9','TNIP1','TPO','UCHL1','UPP1','ZC3H12A','ZEB1','ACPP','MBP','MEFV','NFKB2')

growthExp <- function(model, cell, time)
{
    cycLength <- getCycleLength(model, time, cell)
    return(exp(-1 * cycLength / 48))
}

mitosisExp <- function(model, cell, time)
{
    window <- c(max(time - 2, 0), min(time + 2, model@runTime))

    a1 <- getAxisLength(model, window[1], cell)
    a2 <- getAxisLength(model, window[2], cell)
    if (is.na(a1)) a1 <- 0


    return(ifelse(a2 < a1, 1, 0))
}

SPhaseExp <- function(model, cell, time)
{
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    r1 <- getRadius(model, window[1], cell)
    r2 <- getRadius(model, window[2], cell)

    type <- model@cellTypes[[getCellType(model, time, cell)]]

    return(ifelse(r1 < sqrt(1.5 * type@size) & r2 > sqrt(1.5 * type@size),
        1, 0))
}

contactInhibitionExp <- function(model, cell, time)
{
    return(getLocalDensity(model, time, cell, 3.3))
}

pwyGrowth <- new('Pathway', genes = geneNamesGrowth,
    expressionScale = growthExp)

pwyMitosis <- new('Pathway', genes = geneNamesGtoM,
    expressionScale = mitosisExp, transformMidpoint = 0.05,
    transformSlope = 5 / 0.1)

pwySPhase <- new('Pathway', genes = geneNamesGtoS,
    expressionScale = SPhaseExp, transformMidpoint = 0.05,
    transformSlope = 5 / 0.1)

pwyContactInhibition <- new('Pathway', genes = geneNamesProx,
    expressionScale = contactInhibitionExp)

save(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition,
    file = 'SamplePathways.RData')
