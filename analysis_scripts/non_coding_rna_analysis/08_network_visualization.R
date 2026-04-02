suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(readr)
})

# 1. LOAD DATA
cor_all <- readRDS("ncrna_analysis/All_ncRNA_mRNA_Correlations.rds")
val_targets <- readRDS("ncrna_analysis/Validated_miRNA_mRNA_Targets.rds")

# 2. SELECT TOP INTERACTIONS FOR VISUALIZATION
# We want to show a clear network, so we limit to top hubs
top_mir <- val_targets |>
  arrange(Correlation) |>
  slice_head(n = 50) |>
  select(ncRNA_ID, mRNA_ID, Correlation, biotype)

top_lnc <- cor_all |>
  filter(biotype != "miRNA", abs(Correlation) > 0.7) |>
  arrange(desc(abs(Correlation))) |>
  slice_head(n = 50) |>
  select(ncRNA_ID, mRNA_ID, Correlation, biotype)

edges_df <- bind_rows(top_mir, top_lnc)

# 3. CREATE iGraph OBJECT
net <- graph_from_data_frame(d = edges_df, directed = TRUE)

# 4. SET ATTRIBUTES FOR VISUALIZATION
V(net)$color <- ifelse(V(net)$name %in% edges_df$ncRNA_ID, "gold", "skyblue")
V(net)$size <- ifelse(V(net)$name %in% edges_df$ncRNA_ID, 6, 4)
V(net)$label.cex <- 0.6

# Edge colors based on correlation
E(net)$color <- ifelse(E(net)$Correlation > 0, "firebrick1", "dodgerblue4")
E(net)$width <- abs(E(net)$Correlation) * 2

# 5. SAVE STATIC PLOT
png("ncrna_analysis/ncRNA_mRNA_Network.png", width = 1200, height = 1200, res = 150)
plot(net, 
     layout = layout_with_fr(net),
     vertex.label.dist = 1.0,
     main = "Top ncRNA-mRNA Regulatory Network")
legend("bottomleft", legend=c("ncRNA", "mRNA"),
       col=c("gold", "skyblue"), pch=21, pt.bg=c("gold", "skyblue"),
       pt.cex=1.5, cex=0.8, bty="n", title="Nodes")
legend("bottomright", legend=c("Positive Cor", "Negative Cor (Targeting)"),
       col=c("firebrick1", "dodgerblue4"), lwd=2,
       cex=0.8, bty="n", title="Edges")
dev.off()

# 6. EXPORT NODE/EDGE LISTS FOR EXTERNAL TOOLS (e.g. Cytoscape)
write.csv(as_data_frame(net, what="edges"), "ncrna_analysis/Network_Edges.csv", row.names = FALSE)
write.csv(as_data_frame(net, what="vertices"), "ncrna_analysis/Network_Nodes.csv", row.names = FALSE)

message("Network visualization generated: ncrna_analysis/ncRNA_mRNA_Network.png")
