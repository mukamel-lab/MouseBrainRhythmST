DoHeatmap(ser, features = best_per_clust2, 
          group.bar = T, 
          group.colors = cols_cluster,
          draw.lines = F )+
  guides(col="none")

x=x[,cluster_order]
best_per_clust=apply(x, 2,function(clust1) { which(clust1 == max(clust1))})
best_per_clust2=apply(x, 2,function(clust1) { blah=which(clust1 == max(clust1)); if(blah %in% best_per_clust) { blah=which(clust1[-blah] == max(clust1[-blah]))}; blah})
best_per_clust2=rownames(x)[best_per_clust2]
best_per_clust2=unique(best_per_clust2)


best_per_clust=apply(x, 2,function(clust1) { rev(order(clust1))[1:2]})
best_per_clust2=unique(rownames(x)[as.numeric(best_per_clust)])
