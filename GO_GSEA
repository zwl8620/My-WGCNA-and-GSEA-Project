#在进行富集之前必须搞清楚数据是来自哪个物种，不同物种有不同的基因注释库，本次分析数据来源于Homo sapiens，故使用org.Hs.eg.db这个注释库；
library(org.Hs.eg.db) 
library(clusterProfiler)
library("pathview")

#读取和整理基因ID
gene_symbol=read.table("test0.txt")
#"tur_genes"/"cy_genes"
gene = as.vector(as.matrix(gene_symbol["V1"])) #必须将数据转化为向量或字符串,注意“V1”为含有"SYMBOL"基因名的一列,可用gene_symbol$SYMBOL代替
gene <- unique(gene) #去掉NA值

#GO富集基因ID转换
gene_trans = bitr(gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db") 
#head(gene_trans)

#使用函数enrichGO进行GO富集,"BP"指的是“Biological Process”，“MF”指的是“Molecular Function”，“CC”指的是“Cellular Component”
go_test=enrichGO(as.vector(as.matrix(gene_trans["ENTREZID"])),
                 OrgDb=org.Hs.eg.db,
                 ont="BP",
                 pAdjustMethod="BH",
                 pvalueCutoff=0.05,
                 qvalueCutoff=0.05,
                 minGSSize=1,
                 keyType="ENTREZID",
                 readable=TRUE)
write.table(go_test,file="go_result.xls",sep="\t") #可能有些基因ID会匹配不上

#GO富集结果可视化
#必要时，加上pdf("filename",width = 12, height=9)，因为"BH"名称很长可能让点和柱显示不出来

#绘制有向无环图
pdf("go_test_graph.pdf")
plotGOgraph(go_test,firstSigNodes =10)
dev.off()

#绘制柱状图
pdf("go_test_barplot.pdf")
barplot(go_test, drop=TRUE,showCategory=12)
dev.off()

#绘制气泡图
pdf("go_test_dot.pdf")
dotplot(go_test)
dev.off()
