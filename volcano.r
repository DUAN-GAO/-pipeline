library(EnhancedVolcano)
Args <- commandArgs()
fil <- Args[6]
path <- Args[7]
res <- read.table(fil,
                  sep="\t",
                  header=T,
                  row.names=1,
                  check.names=F,
                  quote="")  
png(path)
keyvals <- rep('green', nrow(res))
names(keyvals) <- rep('Mid', nrow(res))
# fold change > 1.5 & p-value < 0.0001 为高表达
keyvals[which(res$"log2(Fold_change)" > 1.5 & res$"p-value"<0.0001)] <- 'gold'
names(keyvals)[which(res$"log2(Fold_change)" > 1.5 & res$"p-value"<0.0001)] <- 'high'
# fold change < -1.5 & p-value < 0.0001为低表达
keyvals[which(res$"log2(Fold_change)" < -1.5 & res$"p-value"<0.0001)] <- 'royalblue'
names(keyvals)[which(res$"log2(Fold_change)" < -1.5 & res$"p-value"<0.0001)] <- 'low'

EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2(Fold_change)',
    y = 'p-value',
    xlim = c(-4, 4),
    ylim = c(0,15),
    title = '实验组 vs 对照组',
    pCutoff = 10e-3,
    FCcutoff = 1.5,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
    transcriptPointSize = 3.0,
    transcriptLabSize = 3.0,
    colAlpha = 1,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.8,
    colCustom = keyvals,
    border = 'full',
    legend=c('NS','Log (base 2) fold-change','P value',
      'P value & Log (base 2) fold-change'),
    legendPosition = 'right',
    drawConnectors = FALSE,
    boxedlabels = FALSE,
    legendLabSize = 16,
    legendIconSize = 5.0)