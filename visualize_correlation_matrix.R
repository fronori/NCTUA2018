library(tidyverse)

dat <- read_csv("NCTUA2018.csv")
# 2018年センター試験得点の科目間の相関係数
# データ典拠：荘島宏二郎 (2019)
# https://www.ct.u-tokyo.ac.jp/images/c77958b7a6cc887f8f85b225fd31a6f2.pdf

# 倫理や政治・経済と倫理／政治・経済は同時に受験できないようなので分析から除く
df <- dat %>% 
  select(!`倫理／政治・経済`) %>% 
  rename(term = 科目) %>% 
  filter(term != "倫理／政治・経済")

library(corrr)
# tibbleを相関係数行列に変換
mat <- df %>% 
  as_cordf() %>% 
  as_matrix(diagonal = 1)

library(matrixcalc)
# 相関係数行列が対称行列かどうか確認する（データの入力ミスの発見）
mat %>% 
  matrixcalc::is.symmetric.matrix()

# 転置して引き算して全要素が０かを確認してもよい。
diff <- mat - t(mat)
all(diff == 0)

# 相関係数の絶対値を１から引いて距離に変換
dist_one_minus_abs_r <- df %>% 
  mutate(across(!term, ~if_else(is.na(.x), 0, 1 - abs(.x)))) %>% 
  tibble::column_to_rownames("term")
# corrr::as_matrix()を用いても変換できる

# Ward法による階層的クラスタリング
hc <- hclust(d = as.dist(dist_one_minus_abs_r), method = "ward.D2")
# plot(hc)

library(ggdendro)
# ggdendroパッケージの拡張
# #https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
source("Tweaking_ggdendro.R")
hcdata <- dendro_data_k(hc, 4) # とりあえず4個のクラスタに分けてみる

# ラベルを日本語に対応（for Mac)
update_geom_defaults("text", list(family = "HiraKakuPro-W3", colour = "black"))

p1 <- plot_ggdendro(hcdata,
                    direction = "rl",
                    expand.y = 0.2) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.2, 1)) + 
  labs(x = "", y = "Distance = 1 – |r|", 
       title = "階層的クラスタリング: Ward法") + 
  scale_x_reverse() + 
  theme_gray(base_family = "HiraKakuPro-W3")
ggsave("./plot/ggdendro.pdf", plot = p1,
       width = 4, height = 4,
       device = quartz, type = "pdf")

# カラーパレットを指定する準備と確認
library(pals)
par(mfrow = c(5, 1))
pals::pal.bands(pals::brewer.blues)
pals::pal.bands(pals::brewer.bugn)
pals::pal.bands(pals::coolwarm(n = 11)[6:11])
pals::pal.bands(pals::brewer.pubugn)
pals::pal.bands(pals::brewer.ylorrd)
par(mfrow = c(1, 1))

p2 <- df %>% 
  as_cordf(diagonal = 1) %>% 
  stretch() %>% 
  mutate(x = factor(x, levels = hc$labels[hc$order]), 
         y = factor(y, levels = hc$labels[hc$order] %>% rev())) %>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = r)) + 
  scale_fill_gradientn(colours = brewer.blues(11),
                         limits = c(0, 1)) +
  labs(fill = " r", x = "", y = "", 
       title = "相関ヒートマップ") +
  theme_bw(base_family = "HiraKakuPro-W3") +
  theme(aspect.ratio = 1, 
        legend.key.size = unit(0.5, 'cm'), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("./plot/heatmap.pdf", plot = p2,
       width = 5, height = 5,
       device = quartz, type = "pdf")

# corrr::network_plot()のラベルを日本語に対応（for Mac)
library(ggrepel)
update_geom_defaults("text_repel", list(family = "HiraKakuPro-W3"))

# corrr::network_plot()の自己流カスタマイズ
source("customized_network_plot.cor_df.R")

p3 <- df %>% as_cordf() %>% 
  customized_network_plot.cor_df(min_cor = 0.55) + # brewer.bugn, brewer.pubugn
  scale_colour_gradientn(colours = brewer.ylorrd(11), guide = "colourbar", 
                         limits = c(0.5, 1)) + 
  labs(colour = "  r", title = "相関ネットワーク（MDS）") + theme_void(base_family = "HiraKakuPro-W3") + 
  guides(alpha = "none", 
         size = guide_legend(reverse=TRUE, 
                             override.aes=list(colour=brewer.ylorrd(11) %>% .[7:2])))
ggsave("./plot/network_plot.pdf", plot = p3,
  width = 9, height = 7,
  device = quartz, type = "pdf")

library(qgraph)
par(family = "Osaka") # 日本語対応（for Mac）
mat %>% 
  qgraph(layout = "spring", 
         minimum = 0.5, 
         labels = colnames(mat), 
         label.scale = FALSE, 
         node.width = 1.4, 
         title = "相関グラフ（Springレイアウト）")
dev.copy(quartz, file = "./plot/qplot_spring2.pdf", width = 8, height = 8, 
         type = "pdf")
dev.off()




library(glasso)
res_glasso <- glasso(mat, 0.1)
# res_glasso$graphAttributes$Nodes$labels <- res$graphAttributes$Nodes$labels
# res_glasso$graphAttributes$Nodes$names <- res$graphAttributes$Nodes$names
# Plot:
qgraph(res_glasso, layout = "spring")
?glasso
BICgraph <- qgraph(
  CorMat,
  graph = "glasso",
  sampleSize = nrow(bfi),
  tuning = 0,
  layout = "spring",
  title = "BIC",
  details = TRUE
)

adult_zeroorder <- cor(Rogers)

qgraph(adult_zeroorder, layout="spring",
       
       groups = list(Depression = 1:16,
                     
                     "OCD" = 17:26),
       
       color = c("lightblue",
                 
                 "lightsalmon"))