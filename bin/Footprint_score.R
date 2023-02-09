#!/usr/bin/Rscript
Args <- commandArgs()
score=read.table(Args[6],header=F,quote='\t') 
data<-data.frame(pvlaue=sort(c(-score[,5])),sort=1:length(score[,5]),group=c(rep("Non",round(length(score[,5])*0.8)),rep("Yes",length(score[,5])-round(length(score[,1])*0.8))))
text_data <- data.frame( x = c(round(length(data[,1])*0.8),round(length(data[,1])*0.85),round(length(data[,1])*0.9)),y = c(max(data[,1])*0.3,max(data[,1])*0.6,max(data[,1])*0.9),text = c("Top20%", "Top15%","Top10%"))
library("ggplot2")
tiff(Args[7],width = 1400,height = 1400,res=300)
ggplot() +geom_line(data=data, aes(x=data[,2], y=data[,1],color=data[,3]))+geom_point(data=data, aes(x=data[,2], y=data[,1],color=data[,3],fill=data[,3]),size=2, shape=20)+
labs(y="-log10(Footprint score)",x="Footprint")+scale_color_manual(values =c("grey","#9400D3")) +scale_fill_manual(values =c("grey","#9400D3"))+
theme(legend.position = "none",panel.background = element_blank(), panel.border=element_rect(fill='transparent',color="black",size = 1),plot.margin=unit(c(3, 1, 3, 1),"lines"))+
geom_vline(xintercept=c(round(length(data[,1])*0.8),round(length(data[,1])*0.85),round(length(data[,1])*0.9)),col="#778899",linetype="dotted",size=0.8)+
geom_hline(yintercept=c(data[round(length(data[,1])*0.8),1],data[round(length(data[,1])*0.85),1],data[round(length(data[,1])*0.9),1]),col="#778899",linetype="dotted",size=0.8)+
geom_text(data=text_data,aes(x=text_data[,1],y=text_data[,2],label=text_data[,3]),size=2.5,col="#9400D3")
dev.off()
cut_off<-NULL
for (i in seq(1,30,1)) {
cut_off_sub<-data.frame(paste("top",i,"%",sep=""),round(-data[round(length(data[,1])*(1-i*0.01)),1]))
cut_off<-rbind(cut_off,cut_off_sub)
}
colnames(cut_off)<-c("Percentage","score")
write.table(cut_off,Args[8],row.names=FALSE,sep="\t",quote=F)