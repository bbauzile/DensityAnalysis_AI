#=============================================================================
# 			Density analysis paper of palmipeds farms in France
#
#  Ploting incidence and Re curve of avian influenza in South west France
#  created by : Billy Bauzile (billy.bauzile@gmail.com)
#  originally created during PhD: 2019-2022
#
#
#=============================================================================


rm(list=ls())

library("tidyverse")
require(ggpubr)
require(gridExtra)

### incidence #########
workdir = file.path("")

Fday = as.Date("2016-11-28")
dt.epi = read_delim(file.path("Data/data.txt"), delim= " ", col_names = F) |>
  group_by(X7) |>
  filter(X7>-1) |>
  summarise(date2 = X7+Fday,
            n=n()) |>
  ungroup()


lab = c("Reduction \ndensity to the 80th percentile"="F",
        "Reduction \ndensity to the 85th percentile"= "E",
        "Reduction \ndensity to the 90th percentile"="D",
        "Reduction \ndensity to the 95th percentile" ="C",
        "Reduction \ndensity to the 98th percentile" ="B",
        "Baseline"= "A")


rm_farm_scenario= list.files(path =file.path("Data/scenario_density_reduction/"),
                             full.names = T)

no_farm_scenario =lapply(rm_farm_scenario, function(b){
  read_delim(b, delim = " ", col_names = c("INSEE_COM", "palmi_rm", "galli_rm")) %>%
    summarise(scenario= gsub(pattern = ".*/|.txt", "", b) ,
           n_palmi = sum(as.numeric(palmi_rm)),
           scenario_size = 8380 - n_palmi)
})


farm_scenario =  do.call(rbind, no_farm_scenario)

## simulations --incidence
dt.sims =list.files(file.path("Data/simulations/"),pattern = "summari*", full.names = T)

sumfile = lapply(dt.sims, function(d){
  read_csv(d) |>
  mutate(labs = case_when(scenario=="scenario_0_038"~"Baseline",
                          scenario=='scenario_0_01'~"Reduction \ndensity to the 98th percentile",
                          scenario=="scenario_0_0071"~"Reduction \ndensity to the 95th percentile",
                          scenario=="scenario_0_0043"~"Reduction \ndensity to the 90th percentile",
                          scenario=="scenario_0_00345"~"Reduction \ndensity to the 85th percentile",
                          scenario=="scenario_0_00282"~"Reduction \ndensity to the 80th percentile"))

})


incd_data = do.call(rbind,sumfile) |>
  left_join(farm_scenario) |>
  group_by(scenario) |>
  mutate(rel.size = paste0(round(100*as.numeric(str_sub(IC, 1,3)) /scenario_size, 1 ),
                           "%",
                           "(",
                    round(100*as.numeric(str_sub(IC, 6,8))/scenario_size, 1 ),
                    "-",
                    round(100*as.numeric(str_sub(IC, 10, 12))/scenario_size, 1), ")"  ) )


toplot2 = split(incd_data, f= incd_data$scenario)
cols=c( "#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026","#800026")


pl = list()
for( i in seq_along(toplot2)){
  if(isTRUE(i==6)){
    pl[[i]] = toplot2[[6]] %>% as_tibble() %>%
      filter(date_detected>=0, date_detected<117) %>%
      mutate(t= date_detected+ Fday) %>%
      ggplot() +
      geom_ribbon(aes(x = t, ymin = lw_incidence, ymax = up_incidence,
                      fill = labs), alpha = 0.6)+
      geom_ribbon(aes(x = t, ymin = lw_50_incidence, ymax = up_50_incidence,
                      fill = labs))+
      geom_point(data = dt.epi, aes(date2, n),col="#a1d99b", size=2)+
      geom_line(data = dt.epi, aes(date2, n),col="#a1d99b")+
      guides(fill="none", col="none")+
      geom_label(aes(label=paste("Abs.size: ", IC),y=20,x=min(t)+5),
                 hjust = 0,
                 label.size = NA)+
      geom_label(aes(label=paste("Rel.size: ", rel.size),y=15,x=min(t)+5),
                 hjust = 0,
                 label.size = NA)+
      scale_fill_manual("Scenario", values = (cols)[i],
                        aes(labels=labs))+
      guides(fill="none", guide_legend(reverse = T))+
      ggnewscale::new_scale_fill() +
      geom_line(aes(t, incidence ,  col="") , lwd=0.9)+
      scale_color_manual("",values = "black")+
      facet_wrap(labs~., labeller = labeller(labs=lab),
                 strip.position = 'left')+
      scale_x_date("Time",date_labels = "%Y/%m",
                   limits = as.Date(c("2016-11-27", "2017-03-25")))+
      scale_y_continuous("Daily incidence",
                         limits = c(0, 25), breaks = seq(0, 25,10))+
      cowplot::theme_cowplot()+
      theme(text = element_text(face='plain',size = (18/.pt)*3.2),
            strip.background = element_blank(),
            strip.text.y.left= element_text(angle = 0, hjust=1,vjust = 1))
  }
  else{

    pl[[i]] = toplot2[[i]] %>% as_tibble() %>%
      filter(date_detected>-1) %>%
      mutate(t= date_detected+ Fday, date_detected<117) %>%
      ggplot() +
      geom_ribbon(aes(x = t, ymin = lw_incidence, ymax = up_incidence,
                      fill = labs), alpha = 0.6)+
      geom_ribbon(aes(x = t, ymin = lw_50_incidence, ymax = up_50_incidence,
                      fill = labs))+
      geom_label(aes(label=paste("Abs.size: ",IC),y=20,x=min(t)+5),
                 hjust = 0,
                 label.size = NA)+
      geom_label(aes(label=paste("Rel.size: ", rel.size),y=15,x=min(t)+5),
                 hjust = 0,
                 label.size = NA)+
      scale_fill_manual("Scenario", values = (cols)[i],
                        aes(labels=labs))+
      guides(fill="none", guide_legend(reverse = T))+
      ggnewscale::new_scale_color() +
      geom_line(aes(t, incidence),
      linetype= 1,
      col="black",
      lwd=0.9)+
      guides(guide_legend(reverse =T))+
      facet_wrap(labs~., labeller = labeller(labs=lab),
                 strip.position = 'left')+
      scale_x_date("Time",date_labels = "%Y/%m",
                   limits = as.Date(c("2016-11-27", "2017-03-25")))+
      scale_y_continuous("Daily incidence",
                         limits = c(0, 25), breaks = seq(0, 25,10))+
      cowplot::theme_cowplot()+
      theme(text = element_text(face='plain',size = (18/.pt)*3.2,hjust = 0),
            strip.background = element_blank(),
            strip.text.y.left= element_text(angle = 0, hjust=1,vjust = 1)
            )
  }
}

remove_inc_x1 <- cowplot::theme_cowplot(font_size = 14)+
  theme(
    strip.text.y.left= element_text(angle = 0,size = 12, hjust=0,vjust = 1),
  text = element_text(face='plain',size = (14/.pt)*3.2),
  strip.placement = 'outside',
  strip.background = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank()
)
remove_inc_x2 <-cowplot::theme_cowplot(font_size = 14)+
  theme(
  strip.text.y.left= element_text(angle = 0,size = 12, hjust=0,vjust = 1),
  text = element_text(face='plain',size = (14/.pt)*3.2),
  strip.placement = 'outside',
  strip.background = element_blank(),
  # axis.text.x = element_blank(),
  # axis.ticks.x = element_blank(),
  axis.title =  element_blank()
)

sc1= lapply(pl[c(6:2)], function(x) x+remove_inc_x1)
sc2_3 = lapply(pl[c(1)],function(x) x+remove_inc_x2)


p_inc = c(sc1, sc2_3)

inc_p = ggarrange(plotlist = p_inc, ncol = 1, byrow=T)

inc_2 = grid.arrange( grobs=p_inc[1:6], ncol = 1)

save.image(file.path("incidence/incidence_analysis_data.Rdata") )


## re ####
#############

dt.trees = list.files(file.path( "Data/transmission_trees"), "scenario", full.names = T)

dt.trees_doc = lapply(dt.trees, function(ab){
  read_delim(ab) |>
    mutate(scenario= gsub(x=ab, pattern = ".*/|_trees.txt",replacement = ""))
})


type_infected = lapply(dt.trees_doc, function(x){
  relative = x |>
    mutate(farm_type = ifelse(type_infector=="external", "external", "farm")) |>
    group_by(simID, farm_type) |>
    count() |>
    group_by(farm_type) |>
    summarise(Secd = round(median(n),1),
              Secd_lw005 = round(quantile(n, probs = 0.025),1),
              Secd_lw025 = round( quantile(n, probs = 0.25),1),
              Secd_lw075 = round(quantile(n, probs = 0.75),1),
              Secd_lw0975 =round( quantile(n, probs = 0.975),1),
              scenario = unique(x$scenario) ) |>
    summarise(scenario = unique(x$scenario) ,
              farm_type,
              secd = paste0(Secd, "(", Secd_lw025, "-", Secd_lw0975,")" ),
              secd.ratio =paste0(round(Secd/Secd[farm_type=="external"], 1), "(",
                                 round(Secd_lw025/Secd_lw025[farm_type=="external"] ,1),
                                 "-", round(Secd_lw075/Secd_lw075[farm_type=="external"],1) ,")" ) )
})
ratio_secd =  do.call(rbind,type_infected)|>
               
  mutate(labs = case_when(scenario=="scenario_0_038"~"Baseline",
                          scenario=='scenario_0_01'~"Reduction \ndensity to the 98th percentile",
                          scenario=="scenario_0_0071"~"Reduction \ndensity to the 95th percentile",
                          scenario=="scenario_0_0043"~"Reduction \ndensity to the 90th percentile",
                          scenario=="scenario_0_00345"~"Reduction \ndensity to the 85th percentile",
                          scenario=="scenario_0_00282"~"Reduction \ndensity to the 80th percentile")) |>
  select(-scenario) |>
  rename(Scenario = labs) %>%
  filter(!is.na(Scenario),
         farm_type!="external")

dt.RE = read_csv(file.path("RE/re_summarise.csv")) |>
  left_join(ratio_secd, by=c(labs="Scenario"))
#
# cols =  c('#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026', "#800026")
dt.RE$labs = factor(dt.RE$labs,unique(dt.RE$labs)[6:1])
dt.Re =split(dt.RE, dt.RE$labs)

re = list()

for(i in seq_along(dt.Re)){
  re[[i]] = ggplot(dt.Re[[i]])+
  geom_ribbon(aes(x =t, ymin = p_2.5, ymax = p_97.5,
                  fill = labs), alpha = 0.5)+
  geom_ribbon(aes(x = t, ymin = p_25, ymax = p_75,
                  fill = labs))+
    geom_hline(yintercept = 1, lty=2)+
  geom_line(aes(t, p_m),  col="black", lwd=0.6)+
  geom_label(data=dt.Re[[i]] |>
               select(secd, t) |> distinct(),
             aes(label=paste("Secd: ", secd )),
             x=Fday+5,
             hjust = 0,
             y= 3.5, label.size = NA )+
    geom_label(data=dt.Re[[i]] |>
                 select(secd.ratio, t) |> distinct(),
               aes(label=paste("Secd.ratio: ", secd.ratio )),
               x=Fday+5,
               y = 2.7,
               hjust = 0,
               label.size = NA) +
  scale_fill_manual("Scenario", values = (cols)[i])+
  guides(fill="none", guide_legend(reverse = T))+
  scale_x_date("Time",date_labels = "%Y/%m",
               limits = as.Date(c("2016-11-27", "2017-03-25")))+
    scale_y_continuous("", limits = c(0, 4), breaks = seq(0, 4, 2))+
  cowplot::theme_cowplot()+
  theme(text = element_text(face='plain',size = 14,hjust = 0),
        strip.background = element_blank())

}

ggarrange(plotlist = re[6:1], ncol = 1)




# cleaning axis
remove_re_x1 <- cowplot::theme_cowplot(font_size = 14)+
  theme(
  strip.text= element_blank(),
    text = element_text(face='plain',size = (14/.pt)*3.2),
    # strip.placement = 'outside',
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
remove_re_x2 <-cowplot::theme_cowplot(font_size = 14)+
  theme(
    strip.text= element_blank(),
    text = element_text(size = (14/.pt)*3.2),
    # strip.placement = 'outside',
    strip.background = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title =  element_blank()
  )


#
scre= lapply(re[6:2], function(x) x+remove_re_x1)
scre6 =  lapply(re[c(1)],function(x) x+remove_re_x2)

plot_re = c(scre, scre6)
re_p = grid.arrange(grobs =plot_re, ncol = 1)

re_title = expression(R[e])

ggarrange(annotate_figure(inc_2,
                          top = text_grob("Expected daily incidence ", rot = 0,
                                          size=14, hjust = 0.5,
                                           vjust = 0.5),
                          bottom=text_grob("Time",size=14) ),
          annotate_figure(re_p,
                          top = text_grob(re_title, rot = 0,
                                                 size=14, hjust = 0.5,
                                                  vjust = 0.5),
          bottom=text_grob("Time",size=14 ) )
)

ggsave(filename = "RE/combine_plot.tiff",
       dpi = 500, width = 15, height = 12,
        bg= "white",compression = "lzw")




