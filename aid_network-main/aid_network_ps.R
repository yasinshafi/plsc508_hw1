## =============================================================================
## PLSC 508 Problem Set: Foreign Aid Network, 2010
## Data: Oishi, Ito, Murase, Takikawa & Sakamoto (2022), "Evolution of global
##       development cooperation," PLOS ONE 17(8): e0272440.
##       Replication data: https://github.com/kj-oishi/aid_network
##
## Tasks: (1) visualize, (2) centrality, (3) reciprocity + transitivity,
##        (4) community detection.
## Outputs go to output/ (figures as PDF + PNG, tables as CSV, summary.txt).
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
pkgs <- c("igraph", "jsonlite", "ggplot2", "ggrepel", "dplyr")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# igraph 1.x / 2.x compatibility
if (!exists("as_undirected")) as_undirected <- igraph::as.undirected
if (!exists("hits_scores")) hits_scores <- function(graph, weights = NULL)
  list(hub = hub_score(graph, weights = weights)$vector,
       authority = authority_score(graph, weights = weights)$vector)

YEAR  <- 2010          # same cross-section the paper reports in Figs 2-3, 6-7
N_SIM <- 500           # null-model draws for transitivity
set.seed(508)

# dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

## Download replication files if not already present
# base_url <- "https://raw.githubusercontent.com/kj-oishi/aid_network/main/data/"
# for (f in c(sprintf("year_%d_s.json", YEAR), "actors.csv", "LDC.csv")) {
#   dest <- file.path("data", f)
#   if (!file.exists(dest)) download.file(paste0(base_url, f), dest, mode = "wb")
# }

## ---- 1. Build the network ---------------------------------------------------
edges  <- fromJSON(file.path("data", sprintf("year_%d_s.json", YEAR)))
actors <- read.csv("data/actors.csv", stringsAsFactors = FALSE)
ldc    <- read.csv("data/LDC.csv", stringsAsFactors = FALSE)
# LDC.csv holds the full UN LDC list for each year (paper's read_LDC() rule)
ldc_set <- unique(ldc$LDC[ldc$year == YEAR])

# Raw weight = ODA commitment (USD millions). Paper's normalized weight:
# w_ij = f_ij / (s_out_i * s_in_j)  ("gravity" in their code)
edges <- edges |>
  group_by(donor) |> mutate(s_out = sum(weight)) |> ungroup() |>
  group_by(recip) |> mutate(s_in  = sum(weight)) |> ungroup() |>
  mutate(gravity = weight / (s_out * s_in))

g <- graph_from_data_frame(edges[, c("donor", "recip", "weight", "gravity")],
                           directed = TRUE)

V(g)$region   <- actors$region[match(V(g)$name, actors$actor)]
V(g)$region[is.na(V(g)$region)] <- "Unknown"
V(g)$is_donor <- degree(g, mode = "out") > 0
V(g)$is_recip <- degree(g, mode = "in")  > 0
V(g)$role     <- ifelse(V(g)$is_donor & V(g)$is_recip, "Donor & recipient",
                 ifelse(V(g)$is_donor, "Donor", "Recipient"))
V(g)$ldc      <- V(g)$name %in% ldc_set

desc <- data.frame(
  stat  = c("Nodes", "Edges", "Donors (out-degree > 0)",
            "Recipients (in-degree > 0)", "Donor-cum-recipients",
            "Density (directed)", "Total ODA commitments (USD m)"),
  value = c(vcount(g), ecount(g), sum(V(g)$is_donor), sum(V(g)$is_recip),
            sum(V(g)$is_donor & V(g)$is_recip),
            round(edge_density(g), 4), round(sum(E(g)$weight)))
)
write.csv(desc, "output/table0_descriptives.csv", row.names = FALSE)

## ---- 2. Task 1: Visualization ----------------------------------------------
V(g)$strength <- strength(g, mode = "all", weights = E(g)$weight)

reg_levels <- sort(unique(V(g)$region))
pal <- setNames(hcl.colors(length(reg_levels), "Dark 3"), reg_levels)

# Layout: Fruchterman-Reingold on the undirected graph, weighted by the
# paper's normalized "gravity" weight so dependent dyads pull together.
g_lay <- as_undirected(g, mode = "collapse",
                       edge.attr.comb = list(weight = "sum", gravity = "sum"))
lay <- layout_with_fr(g_lay, weights = 1 + 50 * E(g_lay)$gravity / max(E(g_lay)$gravity),
                      niter = 3000)
lay <- apply(lay, 2, function(x) { q <- quantile(x, c(.02, .98))   # pull in outliers
                                   pmin(pmax(x, q[1]), q[2]) })
top_lab <- V(g)$name[order(-V(g)$strength)][1:12]
lab     <- ifelse(V(g)$name %in% c(top_lab, "Bangladesh"), V(g)$name, NA)
lab     <- sub(" \\(People's Republic of\\)", "", lab)

draw_net <- function() {
  par(mar = c(0, 0, 2, 0))
  plot(g, layout = lay,
       vertex.size        = 1.5 + 6 * sqrt(V(g)$strength / max(V(g)$strength)),
       vertex.color       = pal[V(g)$region],
       vertex.shape       = ifelse(V(g)$is_donor, "square", "circle"),
       vertex.frame.color = NA,
       vertex.label       = lab, vertex.label.cex = 0.55,
       vertex.label.color = "black", vertex.label.family = "sans",
       edge.width         = 0.2 + 1.5 * (E(g)$weight / max(E(g)$weight))^0.5,
       edge.color         = adjustcolor("grey40", 0.15),
       edge.arrow.size    = 0)
  title(sprintf("Foreign aid network, %d (ODA commitments)", YEAR), cex.main = 1)
  legend("bottomleft", legend = names(pal), col = pal, pch = 16, cex = 0.6,
         bty = "n", title = "Region / type")
  legend("bottomright", legend = c("Donor", "Recipient"), pch = c(15, 16),
         cex = 0.6, bty = "n", title = "Shape")
}
pdf("output/fig1_network.pdf", width = 10, height = 10); draw_net(); dev.off()
png("output/fig1_network.png", width = 2400, height = 2400, res = 250); draw_net(); dev.off()

## ---- 3. Task 2: Centrality --------------------------------------------------
# Directed, near-bipartite graph -> use direction-appropriate measures.
#  * Donors: out-degree (# recipients) and out-strength (USD).
#  * Recipients: in-degree (# donors) and in-strength (USD).
#  * Hub/authority scores (HITS): donors are hubs, recipients authorities;
#    a recipient scores high if funded by donors who fund many recipients.
hits_w <- hits_scores(g, weights = E(g)$weight)

cent <- data.frame(
  actor      = V(g)$name,
  role       = V(g)$role,
  region     = V(g)$region,
  ldc        = V(g)$ldc,
  out_degree = degree(g, mode = "out"),
  in_degree  = degree(g, mode = "in"),
  out_strength = round(strength(g, mode = "out", weights = E(g)$weight), 1),
  in_strength  = round(strength(g, mode = "in",  weights = E(g)$weight), 1),
  hub        = round(hits_w$hub, 3),
  authority  = round(hits_w$authority, 3)
)

top_donors <- cent |> filter(out_degree > 0) |> arrange(desc(out_strength)) |>
  mutate(rank_strength = row_number()) |>
  arrange(desc(out_degree)) |> mutate(rank_degree = row_number()) |>
  arrange(rank_strength) |> head(10) |>
  select(actor, region, out_strength, rank_strength, out_degree, rank_degree, hub)

top_recips <- cent |> filter(in_degree > 0) |> arrange(desc(in_strength)) |>
  mutate(rank_strength = row_number()) |>
  arrange(desc(in_degree)) |> mutate(rank_degree = row_number()) |>
  arrange(rank_strength) |> head(10) |>
  select(actor, region, ldc, in_strength, rank_strength, in_degree, rank_degree, authority)

write.csv(cent,       "output/centrality_all_nodes.csv", row.names = FALSE)
write.csv(top_donors, "output/table1a_top_donors.csv",   row.names = FALSE)
write.csv(top_recips, "output/table1b_top_recipients.csv", row.names = FALSE)

# Figure: breadth (degree) vs. volume (strength) - are big funders also broad?
cent_plot <- cent |>
  mutate(deg = ifelse(role == "Recipient", in_degree, out_degree),
         str = ifelse(role == "Recipient", in_strength, out_strength),
         side = ifelse(role == "Recipient", "Recipients (in-ties)", "Donors (out-ties)"))
cent_plot <- cent_plot |> group_by(side) |>
  mutate(lab = ifelse(rank(-str) <= 6 | rank(-deg) <= 3 | actor == "Bangladesh",
                      sub(" \\(People's Republic of\\)", "", actor), NA)) |> ungroup()

p2 <- ggplot(cent_plot, aes(deg, str)) +
  geom_point(aes(color = region), alpha = 0.8) +
  geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 30, na.rm = TRUE) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = pal) +
  facet_wrap(~ side, scales = "free_x") +
  labs(x = "Degree (number of partners)", y = "Strength (USD m, log scale)",
       color = "Region / type",
       title = sprintf("Breadth vs. volume of aid ties, %d", YEAR)) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("output/fig2_centrality.pdf", p2, width = 10, height = 6)
ggsave("output/fig2_centrality.png", p2, width = 10, height = 6, dpi = 250)

# Degree vs strength rank correlation (does volume track breadth?)
rho_don <- with(subset(cent, out_degree > 0),
                cor(out_degree, out_strength, method = "spearman"))
rho_rec <- with(subset(cent, in_degree > 0),
                cor(in_degree, in_strength, method = "spearman"))

## ---- 4. Task 3: Reciprocity and transitivity ------------------------------
# 4a. Reciprocity on the directed graph
recip     <- reciprocity(g)
mutual    <- dyad_census(g)
mutual_pairs <- {
  el <- as_edgelist(g)
  key <- paste(el[, 1], el[, 2]); rkey <- paste(el[, 2], el[, 1])
  unique(t(apply(el[key %in% rkey, , drop = FALSE], 1, sort)))
}

# 4b. Transitivity. Donor->recipient ties almost never close triangles
# (a bipartite graph has zero triangles by construction), so ordinary
# transitivity is uninformative here. We report:
#  (i)  Robins-Alexander bipartite clustering: share of 3-paths
#       (d1-r1-d2-r2) closed into 4-cycles, i.e. "if two donors share one
#       recipient, do they share another?" - the bipartite analogue of
#       transitivity;
#  (ii) global transitivity of the donor-donor and recipient-recipient
#       one-mode projections (tie = at least one shared partner), plus
#       Barrat's weighted version (weight = # shared partners).
# Each is compared against a degree-preserving null: directed edge swaps
# (a->b, c->d) => (a->d, c->b) keep every donor's out-degree and every
# recipient's in-degree, so the null holds the bipartite degree sequences fixed.
bip_matrix <- function(graph) {
  A <- as_adjacency_matrix(graph, sparse = TRUE)
  A[Matrix::rowSums(A) > 0, Matrix::colSums(A) > 0]      # donors x recipients
}
project <- function(A) {
  D <- as.matrix(A %*% Matrix::t(A)); R <- as.matrix(Matrix::t(A) %*% A)
  diag(D) <- 0; diag(R) <- 0
  list(donor = graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE),
       recip = graph_from_adjacency_matrix(R, mode = "undirected", weighted = TRUE))
}
ra_clustering <- function(A) {
  S  <- as.matrix(A %*% Matrix::t(A))                   # shared recipients
  c4 <- sum(choose(S[upper.tri(S)], 2))                  # number of 4-cycles
  kd <- Matrix::rowSums(A); kr <- Matrix::colSums(A)
  idx <- Matrix::which(A > 0, arr.ind = TRUE)
  l3 <- sum((kd[idx[, 1]] - 1) * (kr[idx[, 2]] - 1))     # number of 3-paths
  4 * c4 / l3
}
trans_stats <- function(graph) {
  A  <- bip_matrix(graph); pr <- project(A)
  c(bipartite_RA       = ra_clustering(A),
    donor_density      = edge_density(pr$donor),
    donor_trans        = transitivity(pr$donor, type = "global"),
    donor_trans_barrat = mean(transitivity(pr$donor, type = "barrat",
                                           weights = E(pr$donor)$weight), na.rm = TRUE),
    recip_density      = edge_density(pr$recip),
    recip_trans        = transitivity(pr$recip, type = "global"),
    recip_trans_barrat = mean(transitivity(pr$recip, type = "barrat",
                                           weights = E(pr$recip)$weight), na.rm = TRUE))
}
g_bin <- delete_edge_attr(delete_edge_attr(g, "weight"), "gravity")
obs   <- trans_stats(g_bin)

null <- replicate(N_SIM, {
  gr <- rewire(g_bin, with = keeping_degseq(niter = ecount(g_bin) * 10, loops = FALSE))
  trans_stats(simplify(gr))
})

table2 <- data.frame(
  statistic = c("reciprocity (directed)", "transitivity (directed graph)", names(obs)),
  observed  = round(c(recip, transitivity(g, type = "global"), obs), 4),
  null_mean = c(NA, NA, round(rowMeans(null), 4)),
  null_q025 = c(NA, NA, round(apply(null, 1, quantile, 0.025), 4)),
  null_q975 = c(NA, NA, round(apply(null, 1, quantile, 0.975), 4)),
  p_two_sided = c(NA, NA, round(pmin(1, 2 * pmin(
    (rowSums(null >= obs) + 1) / (N_SIM + 1),
    (rowSums(null <= obs) + 1) / (N_SIM + 1))), 3))
)
write.csv(table2, "output/table2_reciprocity_transitivity.csv", row.names = FALSE)

# Figure: observed vs null for the two most informative clustering measures
show <- c(bipartite_RA = "Bipartite (Robins-Alexander) clustering",
          donor_trans  = "Donor-projection transitivity")
null_df <- do.call(rbind, lapply(names(show), function(k)
  data.frame(stat = show[[k]], val = null[k, ])))
obs_df  <- data.frame(stat = unname(show), val = obs[names(show)])
p3 <- ggplot(null_df, aes(val)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(data = obs_df, aes(xintercept = val), color = "firebrick", linewidth = 1) +
  facet_wrap(~ stat, scales = "free") +
  labs(x = "Clustering coefficient",
       y = sprintf("Count (%d degree-preserving rewirings)", N_SIM),
       title = "Observed transitivity (red) vs. degree-preserving null") +
  theme_minimal(base_size = 10)
ggsave("output/fig3_transitivity_null.pdf", p3, width = 9, height = 4)
ggsave("output/fig3_transitivity_null.png", p3, width = 9, height = 4, dpi = 250)

## ---- 5. Task 4: Community detection ---------------------------------------
# Louvain on the undirected donor-recipient graph. Weights = the paper's
# normalized "gravity" weight (rescaled), so communities reflect relative
# dependence (who matters to whom) rather than raw dollar volume, which
# would be dominated by a few large flows.
g_u <- as_undirected(g, mode = "collapse",
                     edge.attr.comb = list(weight = "sum", gravity = "sum"))
w_comm <- E(g_u)$gravity / max(E(g_u)$gravity)

louv_runs <- lapply(1:50, function(i) cluster_louvain(g_u, weights = w_comm))
best      <- louv_runs[[which.max(sapply(louv_runs, modularity))]]
memb      <- membership(best)
# stability: mean NMI between the best partition and the other runs
nmi_stab  <- mean(sapply(louv_runs, function(x) compare(x, best, method = "nmi")))

# robustness: Walktrap and a raw-USD-weighted Louvain
wt       <- cluster_walktrap(g_u, weights = w_comm)
louv_usd <- cluster_louvain(g_u, weights = log1p(E(g_u)$weight))

V(g)$comm <- memb[V(g)$name]

comm_summary <- data.frame(actor = V(g)$name, comm = V(g)$comm,
                           role = V(g)$role, region = V(g)$region,
                           ldc = V(g)$ldc, strength = V(g)$strength) |>
  group_by(comm) |>
  summarise(
    n_nodes      = n(),
    n_donors     = sum(role != "Recipient"),
    n_recipients = sum(role != "Donor"),
    n_ldc        = sum(ldc),
    top_donors   = paste(head(sub(" \\(People's Republic of\\)", "",
                       actor[role != "Recipient"][order(-strength[role != "Recipient"])]), 4),
                       collapse = "; "),
    top_recips   = paste(head(actor[role != "Donor"][order(-strength[role != "Donor"])], 5),
                       collapse = "; "),
    modal_recip_region = names(sort(table(region[role != "Donor"]), decreasing = TRUE))[1],
    .groups = "drop") |>
  arrange(desc(n_nodes))
write.csv(comm_summary, "output/table3_communities.csv", row.names = FALSE)
write.csv(data.frame(actor = V(g)$name, community = V(g)$comm),
          "output/community_membership.csv", row.names = FALSE)

# Region x community cross-tab (recipients only) for the figure
ct <- as.data.frame(table(comm   = V(g)$comm[V(g)$is_recip],
                          region = V(g)$region[V(g)$is_recip]))
p4 <- ggplot(ct, aes(factor(comm), Freq, fill = region)) +
  geom_col() + scale_fill_manual(values = pal) +
  labs(x = "Louvain community", y = "Number of recipients", fill = "Region",
       title = "Regional composition of recipients by community") +
  theme_minimal(base_size = 10)
ggsave("output/fig4b_community_regions.pdf", p4, width = 8, height = 5)
ggsave("output/fig4b_community_regions.png", p4, width = 8, height = 5, dpi = 250)

# Network re-drawn with community colors
comm_pal <- setNames(hcl.colors(max(memb), "Set 2"), seq_len(max(memb)))
draw_comm <- function() {
  par(mar = c(0, 0, 2, 0))
  plot(g, layout = lay,
       vertex.size  = 1.5 + 6 * sqrt(V(g)$strength / max(V(g)$strength)),
       vertex.color = comm_pal[as.character(V(g)$comm)],
       vertex.shape = ifelse(V(g)$is_donor, "square", "circle"),
       vertex.frame.color = NA,
       vertex.label = lab, vertex.label.cex = 0.55, vertex.label.color = "black",
       vertex.label.family = "sans",
       edge.color = adjustcolor("grey40", 0.12), edge.arrow.size = 0,
       edge.width = 0.3)
  title(sprintf("Louvain communities (gravity-weighted), %d; Q = %.3f",
                YEAR, modularity(best)), cex.main = 1)
  legend("bottomleft", legend = paste("Community", names(comm_pal)),
         col = comm_pal, pch = 16, cex = 0.6, bty = "n")
}
pdf("output/fig4a_communities.pdf", width = 10, height = 10); draw_comm(); dev.off()
png("output/fig4a_communities.png", width = 2400, height = 2400, res = 250); draw_comm(); dev.off()

## ---- 6. Summary -------------------------------------------------------------
bgd <- cent[cent$actor == "Bangladesh", ]
sink("output/summary.txt")
cat("=== Foreign aid network,", YEAR, "===\n"); print(desc, row.names = FALSE)
cat("\n--- Centrality ---\n")
cat("Spearman(degree, strength): donors =", round(rho_don, 2),
    "; recipients =", round(rho_rec, 2), "\n")
cat("Bangladesh: in-degree =", bgd$in_degree, "(rank",
    sum(cent$in_degree > bgd$in_degree) + 1, "of", sum(cent$in_degree > 0), ");",
    "in-strength =", bgd$in_strength, "(rank",
    sum(cent$in_strength > bgd$in_strength) + 1, ")\n")
cat("\nTop donors:\n");     print(top_donors, row.names = FALSE)
cat("\nTop recipients:\n"); print(top_recips, row.names = FALSE)
cat("\n--- Reciprocity detail ---\n")
cat("Reciprocity =", round(recip, 4), "; mutual dyads =", mutual$mut, "\n")
if (nrow(mutual_pairs)) print(mutual_pairs)
cat("\n--- Reciprocity / transitivity (observed vs degree-preserving null) ---\n")
print(table2, row.names = FALSE)
cat("\n--- Communities ---\n")
cat("Louvain modularity Q =", round(modularity(best), 3),
    "; k =", max(memb), "; mean NMI across 50 runs =", round(nmi_stab, 3), "\n")
cat("Walktrap Q =", round(modularity(wt), 3), "; k =", length(wt),
    "; NMI with Louvain =", round(compare(wt, best, "nmi"), 3), "\n")
cat("Louvain (log-USD weights) k =", length(louv_usd),
    "; NMI with gravity Louvain =", round(compare(louv_usd, best, "nmi"), 3), "\n\n")
print(as.data.frame(comm_summary), row.names = FALSE)
cat("\nBangladesh is in community", V(g)$comm[V(g)$name == "Bangladesh"], "\n")
sink()

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
message("Done. Outputs in output/")
