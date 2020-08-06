################## Working with "   " as blanks.  -- Current code
m6A_WT_ChIP = read.table("Shen et al. 2016 Table S1 - m6A WT_2.csv", header = T, sep = ",", skip=1, fill = T, quote = "")
m6A_WT_ChIP_Trimmed = m6A_WT_ChIP[(m6A_WT_ChIP$Gene.ID != ""),]
View(m6A_WT_ChIP_Trimmed)

m6A_fip37_ChIP = read.table("Shen et al. 2016 Table S2 - m6A fip37_2.csv", header = T, skip = 1, sep = ",", fill = T, quote = "")
m6A_fip37_ChIP_Trimmed = m6A_fip37_ChIP[m6A_fip37_ChIP$Gene.ID != "",]
View (m6A_fip37_ChIP_Trimmed)

WT_fip37_ChIP_Comparision = merge(m6A_WT_ChIP_Trimmed, m6A_fip37_ChIP_Trimmed, by = intersect("Gene.ID", "Gene.ID"), all = T)
View(WT_fip37_ChIP_Comparision)

#Getting rid of all NA values in Fold.Enrichment in WTchip

ChIP_WTupfip37all = WT_fip37_ChIP_Comparision[!is.na(WT_fip37_ChIP_Comparision$Fold.enrichment..IP.Input..x),]
View(ChIP_WTupfip37all)

# Merging Shen_ChIP to Shen_RNAseq
Shen_RNAseq = read.table("Shen et al. 2016 Table S3 - RNAseq.csv", skip = 1, header = T, sep = ",", fill = T, quote = "")
Z = merge(ChIP_WTupfip37all, Shen_RNAseq, by.x = "Gene.ID", by.y = "Gene.ID", all = T )
View(Z)

Juntawong_TotalDE = read.table("Juntawong Total DE.csv", header = T, sep = ",", skip = 2, fill = T, quote = "")
View(Juntawong_TotalDE)
Juntawong_TotalDE_trimmed
Juntawong_TotalDE_trimmed = Juntawong_TotalDE[!is.na(Juntawong_TotalDE$log2_FC),]
View(Juntawong_TotalDE_trimmed)
#write.table(Juntawong_TotalDE_trimmed, file = "Juntawong_TotalDE_trimmed.csv", quote = F, sep = ",", col.names = NA)

#Only significant up/down in NS/HS Juntawong.
A = Juntawong_TotalDE_trimmed[Juntawong_TotalDE_trimmed$significant != "",]
View(A)

## Merging ChIP_WT/fip37 to TotalDE
WTenrichfip_TotalDE = merge(A, Z, by.x = "AGI", by.y = "Gene.ID", all = T)
View(WTenrichfip_TotalDE)
#write.table(WTenrichfip_TotalDE, file = "WTenrichfrip_TotalDE.csv", quote = F, sep = ",", col.names=NA)

WTenrichfip_TotalDE_2 = WTenrichfip_TotalDE[!is.na(WTenrichfip_TotalDE$significant),]
View(WTenrichfip_TotalDE_2)
WTenrichfip_TotalDE_3 = WTenrichfip_TotalDE_2[!is.na(WTenrichfip_TotalDE_2$Fold.enrichment..IP.Input..x),]
View(WTenrichfip_TotalDE_3)
WT4 = WTenrichfip_TotalDE_3[is.na(WTenrichfip_TotalDE_3$Fold.enrichment..IP.Input..y),]
write.table(WT4, file = "WTenrich_SigDEtotal.csv", quote = F, sep = ",", col.names = NA)
View(WT4)
write.table(WTenrichfip_TotalDE_3, file = "WT3.csv", quote = F, sep = ",", col.names = NA)

######### ChIP Enriched in fip37 NA in WT
X = WT_fip37_ChIP_Comparision[!is.na(WT_fip37_ChIP_Comparision$Fold.enrichment..IP.Input..y),]
View(X)

W = merge(X, Shen_RNAseq, by.x = "Gene.ID", by.y = "Gene.ID", all = T )
View(W)

V = merge(A, W, by.x = "AGI", by.y = "Gene.ID", all = T)
View(V)


U = V[!is.na(V$significant),]
View(WTenrichfip_TotalDE_2)
S = U[!is.na(U$Fold.enrichment..IP.Input..y),]
View(S)
WT4 = WTenrichfip_TotalDE_3[is.na(WTenrichfip_TotalDE_3$Fold.enrichment..IP.Input..y),]
write.table(S, file = "s.csv", quote = F, sep = ",", col.names = NA)
View(WT4)
write.table(WTenrichfip_TotalDE_3, file = "WT3.csv", quote = F, sep = ",", col.names = NA)
