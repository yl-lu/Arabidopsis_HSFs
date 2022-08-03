setwd("C:/project-hsf/RNA-seq_and_ChIP-seq/HSF_TPM.vs.targetGenes_RPKM")

df.hsf.expr <- read.table("meanTPM_HSFs_Col_22C_37C.xls", row.names = 1)
colnames(df.hsf.expr) <- c("TPM_22C30m","TPM_37C30m")
df.hsp.binding <- read.table("HSP70_HSP90_symbol_promoter_2kb_ChIP_RPKM_22C_37C.xls", header = T, row.names = 1)
df.hsf.expr <- log2(df.hsf.expr + 1)
df.hsp.binding <- log2(df.hsp.binding +1)

my_xlim <- c(min(df.hsf.expr),max(df.hsf.expr) + 0.5)
my_ylim <- c(min(df.hsp.binding),max(df.hsp.binding) + 0.2)

head(df.hsf.expr)
head(df.hsp.binding)
hsf.expr.A2.37C <- df.hsf.expr[rownames(df.hsf.expr) %in% c("HSFA2"),]$TPM_37C30m
hsp.binding.A2.37C <- df.hsp.binding[,colnames(df.hsp.binding) %in% c("HSFA2_37C30m")]
df.hsf.expr <- df.hsf.expr[-which(rownames(df.hsf.expr) %in% c("HSFA2")),]
df.hsp.binding <- df.hsp.binding[,-which(colnames(df.hsp.binding) %in% c("HSFA2_37C30m"))]
colnames(df.hsp.binding)
rownames(df.hsp.binding)

pdf("HSF_expression_vs_binding_strength_of_common_HSFs_targets.pdf",
    width = 3*3,
    height = 3*5)
par(mfrow=c(5,3))
for (i in 1:length(rownames(df.hsp.binding))){
  # hsp_22C <- df.hsp.binding$HSFA1a_22C30m[i]
  # hsp_37C <- df.hsp.binding$HSFA1a_37C30m[i]
  df.hsp.binding[,i]
  rownames(df.hsp.binding)[i]
  plot(x = NULL,
       y = NULL,
       xlim = my_xlim,
       ylim = my_ylim,
       type = "n",
       pch = 19,
       lty = 2,
       col = c("steelblue1","firebrick1"),
       xlab = "Expression of HSF log2(TPM+1)",
       ylab = "Promoter ChIP signal",
       main = gsub("_"," ",rownames(df.hsp.binding)[i]))
  legend("topright",
         legend = c("22C30m","37C30m"),
         col = c("steelblue1","firebrick1"),
         pch = 19,
         cex = 0.9,
         bg = NA,
         box.lty = 0)
  points(x = hsf.expr.A2.37C,
         y = hsp.binding.A2.37C[i],
         # y = c(hsp_22C,hsp_37C),
         pch = 19,
         cex = 0.5,
         col = "firebrick1")
  text(x = hsf.expr.A2.37C,
       y = hsp.binding.A2.37C[i] + 0.1,
       labels = "A2",
       # col = "firebrick1",
       srt = +0,
       cex = 0.9)
  for (HSF in rownames(df.hsf.expr)){
    df.hsp.binding.myHSF <- df.hsp.binding[,grep(HSF,colnames(df.hsp.binding))]
    points(x = c(df.hsf.expr[rownames(df.hsf.expr)==HSF,]$TPM_22C30m,
                 df.hsf.expr[rownames(df.hsf.expr)==HSF,]$TPM_37C30m),
           y = c(df.hsp.binding.myHSF[i,1],
                 df.hsp.binding.myHSF[i,2]),
           # y = c(hsp_22C,hsp_37C),
           pch = 19,
           cex = 0.5,
           col = c("steelblue1","firebrick1"))
    lines(x = c(df.hsf.expr[rownames(df.hsf.expr)==HSF,]$TPM_22C30m,
                df.hsf.expr[rownames(df.hsf.expr)==HSF,]$TPM_37C30m),
          y = c(df.hsp.binding.myHSF[i,1],
                df.hsp.binding.myHSF[i,2]),
          lty = 5,
          lwd = 0.1,
          col = "black")
    text(x = df.hsf.expr[rownames(df.hsf.expr)==HSF,]$TPM_37C30m,
         y = df.hsp.binding.myHSF[i,2] + 0.1,
         labels = gsub("HSF","",HSF),
         # col = "firebrick1",
         srt = +0,
         cex = 0.9)
  }
}
dev.off()