#***require packages for the analysis the analysis
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","magick",
"ComplexHeatmap","DECIPHER","philentropy", "tidyr","RCandy")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

list.1<-c("ESI425H","ESI36H","ASM1913465v1","ASM4079078v1","ASM4079076v1","ASM336729v1",
"URI34H","UNY193P","UNY203P","UCT35H","UNY1128P","UNY1083P","UWI248P","UWI247P",
"UCT113H","URI88H","UNY208P","UNY1090P","ASM4079073v1","UWI263P","ASM4079071v1",
"UCT92H","UCT32H","URI47H","UNY149P","UCT110H","URI46H","URI44H","ASM4079079v1",
"URI117H","URI103H","UCT30H","URI86H","URI112H","URI118H","UNY1038P","ASM2466217v1",
"URI33H","URI48H","UNY990P","UWI283P","XYZ459H","UNY172P","UNY1085P","ASM4079075v1",
"B418PP","UCT50H","UNY169P","ASM4079074v1","URI36H","URI56H","B500PP","B331PP","ASM2466215v1",
"PFhe_I_PB_Ill_cons","ESI403H","ESI26H","UCT96H","UCT124H","UNY1032P","URI101H","UCT29H",
"URI87H","URI107H","URI93H","URI102H","UCT109H","URI39H","URI89H","URI41H","URI42H","URI40H",
"URI111H","UCT31H","URI120H","URI91H","ASM4079077v1","ASM2466219v1","ASM4079080v1","ASM215150v1",
"ASM215148v1","ASM215146v1")

list.2<-c("chr","lp5","cp9","lp17","lp21","lp25","cp26","lp28.1","lp28.2","lp28.3",
"lp28.4","cp32.1","cp32.3","cp32.4","cp32.6","cp32.7","cp32.8","cp32.9",
"lp36","lp38","lp54","lp56")

color.1<-c("#16572f","#47b3da","#ef5656")
color.2<-c("#b02b76","#ffa500","#00b2b2","#155832","#4f7b95")
shape.1<-c(15,17,16,19,13)

set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=10, color="black")
	)
# upload and prepare phyloseq objects***
	tab1=read.delim( "recomb_statistics.txt", sep="\t", row.names=1, header=T)
	DT1=data.table(tab1, keep.rownames=T, key="Node")
	
	tab2=read.delim( "sample_data.txt", sep="\t", row.names=1, header=T)
	DT2=data.table(tab2, keep.rownames=T, key="Id")
	
	dend <- ReadDendrogram("SGCs_LD83_rooted.tre")
	
	DT=merge(DT1,DT2, by.x="Node",by.y="Id")
	
	mat1=as.matrix(DT[,c('RST','BAPS_Core_level1','core','accessory','OspC')])
	rownames(mat1)<-DT$Node
	
	mat2=as.matrix(DT[,c('r.m','rho.theta','Bases.in.Recombinations')])
	rownames(mat2)<-DT$Node
	
	DT1=DT[DT$Bases.in.Recombinations != "0" ]
	DT1$value=DT1$Bases.in.Recombinations/DT1$Genome.Length

	tab3=read.delim("tab_recomb_gene.txt", sep="\t", row.names=1, header=T)
	DT3=data.table(tab3, keep.rownames=F, key="CDS")
	
	DT3$Rec.1<-ifelse(DT3$NumRec =="0", "NR","Rec")
	DT3$value.1<-ifelse(DT3$NumRec =="0", "1","1")
	
	
	pp5=ggplot(DT3, aes(x=Repl, y=value.1, fill=Rec.1, color=Rec.1 ))

	pp5$data$Repl <- ordered(pp5$data$Repl, levels=rev(list.2) )

	p5=pp5 + 
	geom_bar(position="fill", stat="identity") +
	scale_y_discrete(limits=c(0,1))+
	theme_bw()+theme_new+ coord_flip()+
	ylab("proportion of recomb.")
	
#	pdf("barplot_recomb_frq.pdf", useDingbats=FALSE)
	print(p5)
#	dev.off()
	DT31=DT3[DT3$NumRec != "0"]
	

	pp6=ggplot(DT31, aes(Repl, NumRec))
	
	pp6$data$Repl <- ordered(pp6$data$Repl, levels=rev(list.2) )
	
	p6=pp6+
	geom_boxplot(outlier.shape=NA)+
	geom_jitter(aes(color=NumAffectedTaxa/82-0.5),size=1.5, width=0.4)+
	scale_color_gradient2(low="#2D5ECD",mid="#E3932C", high="#A3054A")+
	ylab("recomb. events")+ coord_flip()+
	theme_bw()#+theme_new
	
#	pdf("boxplots_recomb_repl.pdf", useDingbats=FALSE)
#	print(p6)
#	dev.off()
	gridExtra::grid.arrange(p5,p6, nrow=2, ncol=2)
	
	
	