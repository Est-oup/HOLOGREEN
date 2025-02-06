#Sample data détaillé 

suivi_palavas <- read.csv(file = "suivi_palavas_R.csv", header = T, sep = ";")
colnames( suivi_palavas) <- c("date", "T.ext", "T.bassins", "pH", "enrich" )
suivi_palavas$date <- suivi_palavas$date - 20000000 
suivi_palavas <- suivi_palavas[suivi_palavas$date <= 210525,] 
# Convertir la colonne "date" en format date
suivi_palavas$date <- ymd(suivi_palavas$date)


#Data utilisées 

envir_data_plot <- data.frame(physeq_16S@sam_data )
envir_data_plot <- envir_data_plot[envir_data_plot$date <= "2021-05-25",]
envir_data_plot$date <- ymd(envir_data_plot$date)


# you'll also want to play with
theme(plot.margin=unit(c(.2,1,.1,1),"cm"))



#Samdata analysis

ga=ggplot_gtable(ggplot_build(ggplot(envir_data_plot, aes(x= date ,color = enrich))+ geom_line(aes(y = ifelse(enrich == "EDM", nitrates, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date > "2021-03-03" & date <= "2021-03-31", nitrates, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-02" & date <= "2021-04-14", nitrates, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-20" & date <= "2021-04-27", nitrates, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-05-04", nitrates, NA_real_)), linetype = "solid")+ 
                                geom_vline(xintercept = c(as.Date("2021-03-03"), as.Date("2021-04-01"), as.Date("2021-04-19"), as.Date("2021-05-03")), color = "indianred2") +
                                geom_point(aes(y=nitrates))+ scale_x_date(date_labels = "%d-%m-%y")+scale_color_manual(values = c("cadetblue3", "indianred2"))+theme_light()+
                                theme(legend.position="top")+ labs(title="", x="", y="[Nitrates]")))#+geom_vline(xintercept=alpha_date_batch)
gb=ggplot_gtable(ggplot_build(ggplot(envir_data_plot, aes(x= date ,color = enrich))+ geom_line(aes(y = ifelse(enrich == "EDM", nitrites, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date > "2021-03-03" & date <= "2021-03-31", nitrites, NA_real_)), linetype = "solhttp://127.0.0.1:47091/graphics/plot_zoom_png?width=1333&height=654id")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-02" & date <= "2021-04-14", nitrites, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-20" & date <= "2021-04-27", nitrites, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-05-04", nitrites, NA_real_), linetype = "solid"), linetype = "solid")+ geom_point(aes(y=nitrites))+ 
                                geom_vline(xintercept = c(as.Date("2021-03-03"), as.Date("2021-04-01"), as.Date("2021-04-19"), as.Date("2021-05-03")), color = "indianred2") +
                                scale_x_date(date_labels = "%d-%m-%y")+scale_color_manual(values = c("cadetblue3", "indianred2"))+theme_light()+
                                theme(legend.position="none")+ labs(title="", x="", y="[Nitrites]")))#+geom_vline(xintercept=alpha_date_batch)
gc=ggplot_gtable(ggplot_build(ggplot(envir_data_plot, aes(x= date ,color = enrich))+ geom_line(aes(y = ifelse(enrich == "EDM", ammonium, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date > "2021-03-03" & date <= "2021-03-31", ammonium, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-02" & date <= "2021-04-14", ammonium, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-20" & date <= "2021-04-27", ammonium, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-05-04", ammonium, NA_real_)), linetype = "solid")+ geom_point(aes(y=ammonium))+
                                geom_vline(xintercept = c(as.Date("2021-03-03"), as.Date("2021-04-01"), as.Date("2021-04-19"), as.Date("2021-05-03")), color = "indianred2") +
                                scale_x_date(date_labels = "%d-%m-%y")+scale_color_manual(values = c("cadetblue3", "indianred2"))+theme_light()+
                                theme(legend.position="none")+ labs(title="", x="", y="[Ammonium]")))#+geom_vline(xintercept=alpha_date_batch)
gd=ggplot_gtable(ggplot_build(ggplot(suivi_palavas, aes(x= date ,color = enrich))+ geom_line(aes(y = ifelse(enrich == "EDM", pH, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date > "2021-03-03" & date <= "2021-03-31", pH, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-02" & date <= "2021-04-14", pH, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-04-20" & date <= "2021-04-27", pH, NA_real_)), linetype = "solid")+
                                geom_line(aes(y = ifelse(enrich == "ENR" & date >= "2021-05-04", pH, NA_real_)), linetype = "solid")+ geom_point(aes(y=pH))+
                                geom_vline(xintercept = c(as.Date("2021-03-03"), as.Date("2021-04-01"), as.Date("2021-04-19"), as.Date("2021-05-03")), color = "indianred2") +
                                scale_x_date(date_labels = "%d-%m-%y")+scale_color_manual(values = c("cadetblue3", "indianred2"))+theme_light()+
                                theme(legend.position="none")+ labs(title="", x="", y="pH")))
ge=ggplot_gtable(ggplot_build(ggplot(suivi_palavas, aes(x= date ))+ geom_line(aes(y=T.bassins))+ geom_point(aes(y=T.bassins))+
                                scale_x_date(date_labels = "%d-%m-%y") + theme_light())) #+geom_vline(xintercept=alpha_date_batch)))

maxWidth = grid::unit.pmax(ga$widths[2:3], gb$widths[2:3], gc$widths[2:3],gd$widths[2:3],ge$widths[2:3])
ga$widths[2:3] <- as.list(maxWidth)
gb$widths[2:3] <- as.list(maxWidth)
gc$widths[2:3] <- as.list(maxWidth)
gd$widths[2:3] <- as.list(maxWidth)
ge$widths[2:3] <- as.list(maxWidth)

grid.newpage()
grid.arrange(arrangeGrob(ga,gb,gc,gd,ge,nrow=5, heights=c(.4,.4,.4,.4,.4)))


grid.arrange(arrangeGrob(ga,gb,gc,nrow=3, heights=c(.4,.3,.3)))


