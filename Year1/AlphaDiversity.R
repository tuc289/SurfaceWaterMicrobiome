# Calculate and visualize alpha diversity.
alpha_S1 <- estimate_richness(Sediment1_otu, measures="invSimpson")
Arb_site1 = sample_data(Sediment1)$Arb_site
alpha_S1 <- cbind(alpha_S1, Arb_site1)

alpha_W1 <- estimate_richness(Water1_otu, measures= "InvSimpson")
Arb_site2 = sample_data(Water1)$Arb_site
alpha_W1 <- cbind(alpha_W1, Arb_site2)

alpha_S2 <- estimate_richness(Sediment2_otu, measures="InvSimpson")
Arb_site3 = sample_data(Sediment2)$Arb_site
alpha_S2 <- cbind(alpha_S2, Arb_site3)

alpha_W2 <- estimate_richness(Water2_otu, measures= "InvSimpson")
Arb_site4 = sample_data(Water2)$Arb_site
alpha_W2 <- cbind(alpha_W2, Arb_site4)

# Test for statistical difference in alpha diversity among samples groups, using Kruskal-Wallis test.
KW_Simpson_S1 <- kruskal.test(InvSimpson ~ sample_data(Sediment1_otu)$Arb_site, data= alpha_S1)
KW_Simpson_S2 <- kruskal.test(InvSimpson ~ sample_data(Sediment2_otu)$Arb_site, data= alpha_S2)
KW_Simpson_W1 <- kruskal.test(InvSimpson ~ sample_data(Water1_otu)$Arb_site, data= alpha_W1)
KW_Simpson_W2 <- kruskal.test(InvSimpson ~ sample_data(Water2_otu)$Arb_site, data= alpha_W2)
  
# Test for statistical differences in alpha diversity among sampling sites, using Kruskal-Wallis multiple comparison test.
install.packages("dunn.test")
library(dunn.test)
dunn_InvSimpson_S1 <- dunn.test(alpha_S1$InvSimpson, sample_data(Sediment1)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_W1 <- dunn.test(alpha_W1$InvSimpson, sample_data(Water1)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_S2 <- dunn.test(alpha_S2$InvSimpson, sample_data(Sediment2)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_W2 <- dunn.test(alpha_W2$InvSimpson, sample_data(Water2)$Arb_site, method="bonferroni",alt=TRUE)
  
# Violin plots (Inverse Simpson by site)
aS1 <- ggplot(alpha_S1, aes(x=Arb_site1, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Red') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS1
aS2 <- ggplot(alpha_W1, aes(x=Arb_site2, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Blue') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS2
aS3 <- ggplot(alpha_S2, aes(x=Arb_site3, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Red') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS3
aS4 <- ggplot(alpha_W2, aes(x=Arb_site4, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Blue') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS4
  
alpha_bac <- subset(alpha_2, amplicon=="16S")
alpha_fun <- subset(alpha_2, amplicon=="ITS")

inv1 <- ggplot(alpha_bac, aes(x=variable, y=value),) + geom_violin(aes(fill=SampleType), trim=FALSE) +
scale_fill_manual(values=c("red","blue")) +
labs(x="Normalization method", y="Inverse Simpson") +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
      text=element_text(size=14, color="black"))+
theme(panel.grid.major = element_blank(), legend.position="none",panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
inv1
inv2 <- ggplot(alpha_fun, aes(x=variable, y=value)) + geom_violin(aes(fill=SampleType), trim=FALSE) +
  scale_fill_manual(values=c("red", "blue"))+
  labs(x="Normalization method", y="Inverse Simpson") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))+
  theme(panel.grid.major = element_blank(), legend.position="none",panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
inv2
  
plot_grid(inv1, inv2, ncol=1, nrow=2, labels=c("A", "B"), label_size=16)

plot_grid(aS1,aS2,aS3,aS4, labels=c("A","B","C","D"), label_size=20)

