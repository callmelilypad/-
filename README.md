# 加载 R 包
library(circlize)
library(RColorBrewer)

# 颜色渐变（红橙色系）
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(15)  # 让颜色更接近

# 正相关的 r 值
r_values_positive <- c(0.999705136, 0.999648582, 0.999073583, 0.997983337, 0.995152304, 
                       0.994631956, 0.994558945, 0.99303204, 0.992743284, 0.991829872, 
                       0.990839141, 0.988467749, 0.981991474, 0.981625132, 0.979852189)

# 基因对名称
gene_pairs_positive <- c("GeneA-PseudoA", "GeneB-PseudoB", "GeneC-PseudoC", 
                         "GeneD-PseudoD", "GeneE-PseudoE", "GeneF-PseudoF", 
                         "GeneG-PseudoG", "GeneH-PseudoH", "GeneI-PseudoI", 
                         "GeneJ-PseudoJ", "GeneK-PseudoK", "GeneL-PseudoL", 
                         "GeneM-PseudoM", "GeneN-PseudoN", "GeneO-PseudoO")

# 拆分母基因和假基因
gene1 <- sapply(strsplit(gene_pairs_positive, "-"), `[`, 1)  # 母基因
gene2 <- sapply(strsplit(gene_pairs_positive, "-"), `[`, 2)  # 假基因

# 创建数据框
df_positive <- data.frame(
  Gene1 = gene1,
  Gene2 = gene2,
  Correlation = r_values_positive
)

# 颜色映射：让颜色都在 OrRd 的范围内
color_mapping <- colorRamp2(range(r_values_positive), c("#FFC300", "#D73027"))  # 黄色到深红色

# **修正 grid.col（必须带名字）**
all_genes <- unique(c(gene1, gene2))  # 提取所有基因
grid.col <- setNames(rep(coul, length.out = length(all_genes)), all_genes)  # 让颜色和基因一一对应

# **调整圆的大小**
circos.clear()
circos.par(
  track.margin = c(0.01, 0.01),  # **减少扇形块之间的间距**
  canvas.xlim = c(-1.2, 1.2),    # **让整个圆缩小一点**
  canvas.ylim = c(-1.2, 1.2)
)

# **绘制和弦图**
chordDiagram(
  x = df_positive[, 1:2],  # 只取 Gene1 和 Gene2 列
  grid.col = grid.col,      # **每个基因都有自己的颜色**
  transparency = 0.3,       # 透明度适中
  col = color_mapping(df_positive$Correlation),  # **连线颜色按 r 值**
  annotationTrack = "grid"  # **显示基因名称**
)

# **添加基因名称**
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name = get.cell.meta.data("sector.index")  # 获取当前基因名称
  xlim = get.cell.meta.data("xlim")  # 获取当前扇形块范围
  ylim = get.cell.meta.data("ylim")  # 获取当前扇形块Y轴范围
  
  # **在扇形块中部显示名称**
  circos.text(mean(xlim), ylim[1] + 0.5, sector_name, 
              facing = "clockwise",  # 让文字沿着圆周排列
              niceFacing = TRUE,      # 让文字方向更自然
              adj = c(0, 0.5),        # 调整对齐方式
              cex = 0.6)              # **缩小字体**
})

# 添加标题
title("Chord Diagram of Positively Correlated Gene-Pseudogene Pairs") 
# 在右侧添加竖直渐变色图例（类似热图图注）
library(ComplexHeatmap)  # 需要加载此包

lgd <- Legend(
    col_fun = color_mapping, 
    title = "r-value",
    at = c(0.98, 0.99, 1.0),  # 设置要显示的刻度值（根据实际范围调整）
    labels = c("0.98", "0.99", "1.00"),  # 自定义标签显示
    direction = "vertical",    # 竖直方向
    legend_height = unit(4, "cm")  # 调整图例高度
)

# 将图例绘制在右侧空白处
draw(lgd, x = unit(0.92, "npc"), y = unit(0.7, "npc"))  # 调整x/y控制位置
