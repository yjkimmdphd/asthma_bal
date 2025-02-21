
###
# 7. (Alternative) Module merging
### 


ME.dissimilarity = 1-cor(MEs, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic

png(file.path(output_folder,"METree.png"), width = 800, height = 600)

plot(METree)
abline(h=0.18, col = "red")
abline(h=0.25, col = "red")  #a height of h corresponds to correlation of 1.00 - h (i.e.,  all of the modules which are more than 85% similar if h=0.15.)
dev.off()

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = 0.18) # merge the modules which are below the threshold

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# plot dendrogram showing both the orginal and merged module colors

png(file.path(output_folder,"Gene_dendrogram_and_module_colors_for_original_and_merged_modules_h-018.png"), width = 800, height = 600)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

write.table(mergedColors,file.path(output_folder,"mergedColors_h-018.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)
write.table(mergedMEs,file.path(output_folder,"mergedMEs_h-018.txt"), sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)