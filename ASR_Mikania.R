#The script below is freely available and does not offer any guarantee. The use of it is entirely the responsibility of the
#user, as well as any possible eventuality resulting from its misuse.
#This code was tested using R version 3.6.3 and castor 1.6.6, on Windows 10.

#The script was created for phylogenetic analyzes using the "castor" package in order to maintain only the essentials for
#its execution and visualization of the results. Several notes were made in order to help the user. Curly braces "{}" were
#used when necessary to speed up the process. A loop in parallel computing was applied in some cases to automate certain
#processes, however, individual execution is also available. The script was built so that the user can make as few manual
#adjustments as possible, however, feel free to improve it and share it with the community.

#Script author: Fábio Alves.

{memory.limit(size = 2*memory.limit()) #Increased limit of available memory
Packages <- c("ape", "dplyr", "mgsub", "castor", "strap", "phyloch", "phytools", "geoscale", "doParallel")
lapply(Packages, library, character.only = T)
numCores <- detectCores() #Detect the number of CPU cores on the current host. It is important for parallel computing
registerDoParallel(numCores)} #Register the parallel backend with the foreach package

#Read the tree for analysis ("tree" object)
#Read the caracter states document ("tip.states" and "tip_states" objects. Both have different purposes along the script)
{tree = read_tree( file = "./mikania.tre",
                 edge_order = "cladewise",
                 include_edge_lengths = T,
                 look_for_edge_labels = T,
                 look_for_edge_numbers = T,
                 include_node_labels = T,
                 underscores_as_blanks = F,
                 check_label_uniqueness = F,
                 interpret_quotes = F,
                 trim_white = T)
  
tip.states = read.table("./mikania_states.csv", h = T, sep = ";", dec = ",")
tip_states = lapply(tip.states, as.numeric)
Ntips = length(tree$tip.label) #Number of tips
Nnodes = tree$Nnode #Number of nodes
root_age = get_tree_span(tree)$max_distance} #Root age
#node.depth(tree) #Return the depths of nodes and tips. Count number of species of nodes and tips
#branching.times(tree) #Computes the branching times of a tree, that is the distance from each node to the tips


#Plot tree with posterior probabilities, ages 95% HPD and geologic time scale
tree.plotter <- function(basepath, xmin, xmax) {
  base_tree <- read.nexus(paste(basepath, "tree_beast.tre", sep = ""))
  base_tree$root.time <- max(nodeHeights(base_tree))
  annot_tree <- read.beast(paste(basepath, "tree_beast.tre", sep = ""))
  annot_tree$root.time <- max(annot_tree$height)
  age_table <- read.table(paste(basepath, "age_ranges.txt", sep = ""), stringsAsFactors = F)
  params <- read.table(paste(basepath, "loganalyser_params.txt", sep = ""), header = T, stringsAsFactors = F)
  if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
    annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
    annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
  } else {
    annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
    annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
  }
  if (length(params$mean[params$statistic == "offset"] != 0)) {
    offset <- params$mean[params$statistic == "offset"]
    base_tree$root.time <- base_tree$root.time + offset
    annot_tree$min_ages <- annot_tree$min_ages + offset
    annot_tree$max_ages <- annot_tree$max_ages + offset      
  } else {
    base_tree$root.time <- base_tree$root.time
    annot_tree$min_ages <- annot_tree$min_ages
    annot_tree$max_ages <- annot_tree$max_ages
  }
  age_mat <- cbind(as.numeric(age_table[,2]), as.numeric(age_table[,3]))
  rownames(age_mat) <- age_table[,1]
  colnames(age_mat) <- c("FAD", "LAD")
  geoscalePhylo(ladderize(base_tree, right = T), ages = age_mat, show.tip.label = F, x.lim = c(xmin, xmax), y.lim = c(2, 77),
                units = c("Period", "Epoch", "Age"), cex.tip = 0.5, cex.age = 1, cex.ts = 1, width = 1.5, tick.scale = 5)
  nodelabels(round(annot_tree$posterior, 2), adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "firebrick2", font = 2)
  tiplabels(sub("_", " ", base_tree$tip.label), adj = 0, frame = "n", cex = 0.6, col = "black", font = 3, offset = 0.1)
  T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  for(i in (Ntip(base_tree) + 1):(base_tree$Nnode + Ntip(base_tree))) {
    lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(base_tree)],
                T1$root.time - annot_tree$max_ages[i - Ntip(base_tree)]),
          y = rep(T1$yy[i], 2), lwd = 6, lend = 0,
          col = make.transparent("blue", 0.4))}}

tree.plotter("./Mikania/", -0.8, 11)

#Export tree plot
tiff("TREE.tiff", units="px", width=6000, height=3375, res=300)
tree.plotter("./Mikania/", -0.8, 11)
dev.off()


#Phylogenetic Autocorrelation Function Of A Numeric (continuous) Trait ####
ACF = get_trait_acf(tree, tip_states$chromosome, Npairs=1e7, Nbins=6)
plot(ACF$distances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")


#Trait Depth (All) - Phylogenetic signal ####
#Parallel loop process for calculating the trait depth of each character state
trait_depth_all = foreach(i = seq_along(tip_states), .packages = "castor", .combine = cbind) %dopar% {
trait_depth = get_trait_depth(tree,
                              tip_states[[i]],
                              min_fraction = 0.9,
                              count_singletons = T,
                              singleton_resolution= 0,
                              weighted = F,
                              Npermutations = 1000000)}

#Write .csv document with all trait depths
write.table(trait_depth_all, "./Mikania/trait_depth.csv", sep = ";", dec = ".")

#Trait Depth (Individual) - Phylogenetic signal ####
#Process for calculating the trait depth of one character state at a time
trait_depth = get_trait_depth(tree,
                              tip_states$Aragorn, #Select here the column corresponding to the desired character state
                              min_fraction = 0.9,
                              count_singletons = T,
                              singleton_resolution= 0,
                              weighted = F,
                              Npermutations = 1000000)

#Print the mean depth and standard deviation of the individual trait depth
#cat(sprintf("Mean depth = %g, std = %g\n",trait_depth$mean_depth,sqrt(trait_depth$var_depth)))

#Plot Trait Depth ####
#$positive_clades #Indices of tips and nodes (from 1 to Ntips+Nnodes) that were found to be positive in the trait
#You need to change the title and some limits according to your tree tip/node value. The rest is optional.
#The green number shows the positive species for the corresponding node. The number in red is just a tip/node ID
{plot.phylo(ladderize(tree), show.tip.label = F)
title(main = "Hab.")
legend("left", legend = c("Mean Depth", trait_depth$mean_depth, "P Value", trait_depth$P), box.col = "white")
nodelabels(node = trait_depth$positive_clades, adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
tiplabels(text = tree$node.label[trait_depth$positive_clades], trait_depth$positive_clades,
          pch = 19, frame = "n", cex = 0.7, col = "black")
tiplabels(text = trait_depth$positives_per_clade[trait_depth$positive_clades], trait_depth$positive_clades, 
          adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "darkgreen", font = 2)
b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
    filter(as.data.frame(trait_depth$positives_per_clade) == 1 %in% 
    row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. < 133)))
tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
          adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
    filter(as.data.frame(trait_depth$positives_per_clade) == 0 %in% 
    row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. < 133)))
tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
          adj = 0, frame = "n", cex = 0.7, col = "black", font = 3, offset = 0.2)
axisPhylo(1)}

#$positives_per_clade #Number of descending tips per clade (tip or node) that were positive in the trait
#Use "!" before the filter, inverts the selection
{plot.phylo(ladderize(tree), tip.color = "royalblue1", cex = 0.7, label.offset = 0.2, y.lim = 77, x.lim = 12)
t = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(trait_depth$positives_per_clade)))) %>% 
    filter(as.data.frame(trait_depth$positives_per_clade) > 0 %in% 
    row.names(as.data.frame(trait_depth$positives_per_clade))) %>% filter(. > 79)))

nodelabels(node = t, adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
tiplabels(text = tree$node.label[t], t, pch = 1, frame = "n", cex = 0.7, col = "black")
tiplabels(text = trait_depth$positives_per_clade[t], t, 
          adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "darkgreen", font = 2)
axisPhylo(1)}


#Ancestral State Reconstruction (ASR) step ####
#Map states of a discrete trait to integers (All)
#Parallel loop process for map all states
state_map_all = foreach(i = seq_along(tip.states), .packages = "castor", .combine = rbind) %dopar% {
  state_map = map_to_state_space(tip.states[[i]])}

#ASR with Markov (Mk) models (All)
#Parallel loop process for calculating the ASR of each character with all known tip states
amk_all = foreach(i = c(13, 20, 21), .packages = "castor", .combine = rbind) %dopar% {
  amk = asr_mk_model(tree, state_map_all[,3][[i]], state_map_all[,1][[i]], rate_model = "ARD",
                     include_ancestral_likelihoods = T, store_exponentials = T, Ntrials = numCores, Nthreads = numCores)}

#ASR with hidden state prediction via Mk model max-likelihood (All)
#Parallel loop process for calculating the ASR of each character with tip states absence, i.e. hidden state
hmk_all = foreach(i = c(2, 5, 8, 14), .packages = "castor", .combine = rbind) %dopar% {
  hmk = hsp_mk_model(tree, state_map_all[,3][[i]], state_map_all[,1][[i]], rate_model = "ARD", include_likelihoods = T,
                     store_exponentials = T, Ntrials = numCores, Nthreads = numCores)}


#Map states of a discrete trait to integers (Individual)
#Process for map one state at a time. You need to specify the character column from tip.states object changing the number
state_map = map_to_state_space(tip.states[,2])

#ASR with Mk models and rerooting (Individual)
#Process for calculating the ASR of one character at a time with all known tip states
{amk = asr_mk_model(tree, state_map$mapped_states, state_map$Nstates, rate_model = "ARD",include_ancestral_likelihoods = T,
                    store_exponentials = T, Ntrials = numCores, Nthreads = numCores)$ancestral_likelihoods
amk_estimated = max.col(amk[1:(Nnodes),])
amk_estimated} #Print estimated node states

#ASR with hidden state prediction via Mk model max-likelihood (Individual)
#Process for calculating the ASR of one character at a time with tip states absence, i.e. hidden state
{hmk = hsp_mk_model(tree, state_map$mapped_states -1, state_map$Nstates -1, rate_model = "ARD", include_likelihoods = T,
                   store_exponentials = T, Ntrials = numCores, Nthreads = numCores)$likelihoods
hmk_estimated = max.col(hmk[1:(Ntips + Nnodes),])
hmk_estimated} #Print estimated tip and node states


#Plot ASR for AMK approach with pie chart likelihoods and geologic time scale
amk.plotter <- function(basepath, xmin, xmax) {
plot.phylo(ladderize(tree), show.tip.label = F)
colors = c("royalblue1", "firebrick2", "gold", "darkgreen", "sienna")
nodelabels(round(amk, 2), adj = c(0.5, 0.5), frame = "n", cex = 0.25, col = "#00000000", font = 2, piecol = colors, pie = amk)
tiplabels(sub("_", " ", tree$tip.label), adj = 0, frame = "n", cex = 0.6, col = "black", font = 3, offset = 0.001)
rr = data.frame(V1 = map_to_state_space(tip.states[,12])$mapped_states -1,
                V2 = map_to_state_space(tip.states[,11])$mapped_states -1)
                #V3 = map_to_state_space(tip.states[,7])$mapped_states -1)
tiplabels(rr, adj = c(0.5, 0.5), frame = "n", cex = 0.12, col = "#00000000", font = 2, piecol = colors, pie = rr)
legend(x = "bottomleft", title = "Deeply lobed leaf blade", legend = state_map$state_names, pch = 21, pt.cex = 2.5, pt.bg = colors,
       cex = 1.5, bty = "n")}

amk.plotter("./Mikania/", -0.8, 11)

#Export AMK plot
tiff("AMK.tiff", units = "px", width = 6000, height = 3375, res = 300)
amk.plotter("./Mikania/", -0.8, 11)
dev.off()


#Plot ASR for HMK approach with pie chart likelihoods and geologic time scale
hmk.plotter <- function(basepath, xmin, xmax) {
plot.phylo(ladderize(tree), show.tip.label = F)
colors = c("royalblue1", "firebrick2", "gold", "darkgreen", "sienna")
nodelabels(round(hmk[((max(tree$Nnode)+20):max(tree$edge)),1], 2), adj = c(0.5, 0.5), frame = "n", cex = 0.25,
           col = "#00000000", font = 2, piecol = colors, pie = hmk[((max(tree$Nnode)+20):max(tree$edge)),])
tiplabels(sub("_", " ", tree$tip.label), adj = 0, frame = "n", cex = 0.6, col = "black", font = 3, offset = 0.001)
tiplabels(round(hmk[(1:(max(tree$Nnode)+19)),1], 2), adj = c(0.5, 0.5), frame = "n", col = "#00000000", font = 3,
          cex = 0.12, piecol = colors, pie = hmk[(min(tree$edge):(max(tree$Nnode)+19)),])
legend(x = "bottomleft", title = "Cypsela", legend = state_map$state_names[-1], pch = 21, pt.cex = 2.5, pt.bg = colors,
       cex = 1.5, bty = "n")}

hmk.plotter("./Mikania/", -0.8, 11)

#Export HMK plot
tiff("HMK.tiff", units = "px", width = 6000, height = 3375, res = 300)
hmk.plotter("./Mikania/", -0.8, 11)
dev.off()


#Plot ASR for AMK approach
#You need to change the title, legend and some limits according to your tree tip/node value. The rest is optional.
{plot.phylo(ladderize(tree), y.lim = 77, x.lim = 12, show.tip.label = F)
title(main = "A.S.R.: Habit")
legend("left", legend = c("C = Caulirosula", "S = Shrub", "T = Trees"),
       text.col = c("firebrick2", "royalblue1", "darkgreen"), box.col = "white")
#nodelabels(((Ntips + 1):(Ntips + Nnodes)), adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
tiplabels(tip.states[,12], adj = 0, frame = "n", cex = 0.6, col = "darkgreen", font = 2)
nodelabels(mgsub(amk_estimated[1:Nnodes], state_map$name2index, state_map$state_names),
           adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "purple3", font = 2)
b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(state_map$mapped_states)))) %>% 
    filter(as.data.frame(state_map$mapped_states) == 1 %in% 
    row.names(as.data.frame(state_map$mapped_states)))))
tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
          adj = 0, frame = "n", cex = 0.7, col = "firebrick2", font = 3, offset = 0.2)
c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(state_map$mapped_states)))) %>% 
    filter(as.data.frame(state_map$mapped_states) != 1 %in% 
    row.names(as.data.frame(state_map$mapped_states)))))
tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
          adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
axisPhylo(1)}


#Plot ASR for HMK approach
#You need to change the title, legend and some limits according to your tree tip/node value. The rest is optional.
{plot.phylo(ladderize(tree), y.lim = 77, x.lim = 12, show.tip.label = F)
title(main = "A.S.R.: Pappus type")
legend("left", legend = c("S = Setose", "P = Paleaceous"), text.col = c("firebrick2", "royalblue1"), box.col = "white")
#nodelabels(((Ntips + 1):(Ntips + Nnodes)), adj = c(0.5, -0.5), frame = "n", cex = 0.6, col = "red")
tiplabels(mgsub(hmk_estimated[1:Ntips], state_map$name2index -1, state_map$state_names),
          adj = 0, frame = "n", cex = 0.6, col = "darkgreen", font = 2)
nodelabels(mgsub(hmk_estimated[80:(Ntips + Nnodes)], state_map$name2index -1, state_map$state_names),
           adj = c(1.1, 1.1), frame = "n", cex = 0.6, col = "purple3", font = 2)
b = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(hmk_estimated)))) %>% 
    filter(as.data.frame(hmk_estimated) == 1 %in% 
    row.names(as.data.frame(hmk_estimated))) %>% filter(. < 80)))
tiplabels(text = sub("_", " ", tree$tip.label[b]), b, 
          adj = 0, frame = "n", cex = 0.7, col = "royalblue1", font = 3, offset = 0.2)
c = as.numeric(as.matrix(as.data.frame(as.numeric(row.names(as.data.frame(hmk_estimated)))) %>% 
    filter(as.data.frame(hmk_estimated) != 1 %in% 
    row.names(as.data.frame(hmk_estimated))) %>% filter(. < 80)))
tiplabels(text = sub("_", " ", tree$tip.label[c]), c, 
          adj = 0, frame = "n", cex = 0.7, col = "firebrick2", font = 3, offset = 0.2)
axisPhylo(1)}
