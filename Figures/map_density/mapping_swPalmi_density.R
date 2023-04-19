#=============================================================================
# 			Density analysis paper of palmipeds farms in France
#
#  Mapping density of palmipeds farms in South west France
#  created by : Billy Bauzile (billy.bauzile@gmail.com)
#  originally created during my PhD: 2019-2022
#
#
#=============================================================================

library(sf)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(magrittr)

workdir = file.path("Density_paper/map_density")





#department of the sw..selected based on the department that were in 100km vicinity of
# farms that were detected during 2016-2017
dt_map_dp = st_read(file.path(workdir, "../Data/Map_data/shapefile_national_mapFR.shp"))


#southwest department
sw = c("33", "40","64","65","32", "47", "24", "46", "82","81","24", "31",
       "11", "12", "34","48", "09",  "16","17", "19", "15", "30")

#communes
dt_map_cm =  st_read(file.path(workdir, "../Data/Map_data/shapefile_southwFR_map.shp") )


dep_nom = dt_map_dp %>%
  dplyr::filter(CODE_DEPT %in% sw) %>%
  dplyr::select(NOM_DEPT, CODE_DEPT) %>%
  arrange(CODE_DEPT)%>%
  st_drop_geometry() %>%
  as_tibble() %>%
  mutate(study_area = paste0(str_to_title(NOM_DEPT), ",",CODE_DEPT))

cnames <- c("id", "department", "insee", "type", "x", "y",
            "date_detection", "date_culling", "symptomatic",
            "preventively_culled", "delay_empty",
            "op", "is_pag", "in_zrp")


#epidata from the model
epi16 <- read_delim(file.path(workdir, "../Data/data.txt"), delim = " ", col_names = cnames)

epi_palmi = epi16 %>%
  dplyr::select(id,insee, type) %>%
  dplyr::filter(type=="palmipedes") %>%
  group_by(insee) %>%
  count() %>%
  rename(n_baseline =n)

# palmipeds farms density in km2
# This first file created by counting the number of duck farms in each communes and divide it by the superficie of the commune.
# Since surperfie in the shape file is hectare, we convert it to km2
palmi_dens = read_csv(file.path(workdir, "../Data/dt_palmi_density_sw.csv")) %>%
  dplyr::mutate(duck_dn2= duck_dn*100,
                lab = cut(duck_dn2,
                          breaks = c( 0.282, 0.345, 0.43, 0.71, 1.0, 3.80)),
                INSEE_COM = ifelse(str_count(INSEE_COM)<5,
                                   paste0("0", INSEE_COM), INSEE_COM) ) %>%
  dplyr::select(dpt, INSEE_COM, duck_dn, duck_dn2, lab)


##base plot
mp_sw2 = merge(dt_map_cm, epi16 %>%
                 dplyr::filter(date_detection>-1) %>%
                 dplyr::select(id, type, insee),
               by.x="INSEE_COM", by.y="insee", all=T)

mp_sw2v = st_as_sf(mp_sw2)


geom= st_geometry(mp_sw2[!is.na(mp_sw2$id), ] )

x=st_union(st_buffer(geom, 100000), by_feature = F) #100km
xx= st_intersection(x,  dt_map_dp, sparse = T) #epi area

#
dt_map_dp$epi="1"

tts= ggplot()+
  geom_sf(data=dt_map_dp %>%
            dplyr::filter(CODE_DEPT %in% sw), aes(fill=factor(epi)  )) +
  geom_sf(data=xx, fill="grey95") +
  scale_fill_manual(values = c("#2B8CBE","grey95"),
                    name = "", #"Epizootic area",
                    labels = c("Outside\nepidemic area"),
                    guide = guide_legend(aes = list(alpha = 0.4),
                                         order = 1,
                                         nrow=2, byrow = T, size= 14))+
  theme_map()+
  #add space to the right of the legend
  theme(legend.text = element_text(size=16,
                                   margin = margin(r = 3, unit = "cm") ),
        legend.position = "top",
        plot.margin = unit(c(1, 1, 0,0), "cm")
  )

tts
## Preparing the density data
# col scale
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
cols <- f("YlOrRd") #[rev(c(9,7,5,3,2, 1))])

cols2 = (c('#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026', "#800026") ) ##initial palette

# cols2 =cols
nom =  expression(atop(Commune-level~density~of~duck~farms,(farms/km^{2}) ))


bn= c(2.82, 3.45,4.3, 7.1, 10, 38.0)/10

############

mp_sw_dn = dt_map_cm %>%
  dplyr::select(INSEE_COM, SUPERFICIE) %>%
  st_drop_geometry() %>%
  mutate(SUPERFICIE_km= as.numeric(SUPERFICIE)/100) #superficie in km2

## adding the density to the map data
test2 =lapply(bn, function(a){
  # a=bn[6]
  tmp3 = epi_palmi  %>%
    dplyr::inner_join(palmi_dens %>%
                        dplyr::select(INSEE_COM, duck_dn) %>%
                        dplyr::filter(!is.na(duck_dn)),
                      by=c("insee"="INSEE_COM")) %>%
    mutate(scenario_palmi=trunc(case_when(duck_dn*100>= a~
                                            n_baseline * ( a/(duck_dn*100) ) ))) %>%
    mutate(scenario_palmi=ifelse(is.na(scenario_palmi), n_baseline, scenario_palmi)) %>%
    mutate(scenario = ifelse(isTRUE(a==3.800), "baseline", a))%>%
    group_by(scenario) %>%
    left_join(mp_sw_dn, by=c("insee"="INSEE_COM")) %>%
    mutate(palmi2 = (scenario_palmi/SUPERFICIE_km) ) %>%
    mutate(palmi2 = ifelse(palmi2 >a, a, palmi2 ) )

  tmp3$lab = cut(tmp3$palmi2, c(0, 0.282, 0.345, 0.43, 0.71, 1, 3.80),
                 include.lowest = F, ordered_result = T)

  tmp4 = tmp3  %>%
    mutate(td_lab = case_when(scenario=="baseline"~"Baseline",
                              scenario=="1.000"~"Reduction \ndensity to the 98th percentile",
                              scenario=="0.710"~"Reduction \ndensity to the 95th percentile",
                              scenario=="0.430"~"Reduction \ndensity to the 90th percentile",
                              scenario=="0.345"~"Reduction \ndensity to the 85th percentile",
                              scenario=="0.282"~"Reduction \ndensity to the 80th percentile"))

  return(tmp4)
})

#linking the scenario with the map data
dt.dens = lapply(test2, function(x){
  x %>%
    dplyr::select(insee,n_baseline,scenario_palmi, duck_dn, td_lab, palmi2) %>%
    ungroup() %>%
    distinct(insee,.keep_all=T)
})

scenario_den = do.call(rbind, dt.dens)

# write_csv(scenario_den, "density_scenario.csv")


dt.map = lapply(test2, function(yy){
  # yy= test2[[6]]
  tmp5 = yy %>%
    merge(dt_map_cm, by.x= "insee", by.y="INSEE_COM", all=T) %>%
    st_as_sf()

  tmp6 = st_intersection(tmp5, xx)

  # tmp6$lab2 = (tmp6$lab)
  # tmp6$lab2[is.na(levels(tmp6$lab2) )]= "0"
  tmp6 = tmp6 %>%
    mutate(col_lab2 = case_when(lab == "(0,0.282]"~"(0,0.282]",
                                lab == "(0.282,0.345]"~"(0.282,0.345]",
                                lab == "(0.345,0.43]"~"(0.345,0.43]",
                                lab == "(0.43,0.71]"~"(0.43,0.71]",
                                lab == "(0.71,1]"~"(0.71,1]",
                                lab == "(1,3.8]"~"(1,3.8]",
                                is.na(lab)~"0"))
  niveau =  str_sort(tmp6$col_lab2) %>% unique()
  niveau = niveau[c(length(niveau), 1: (length(niveau)-1) )]
  tmp6$col_lab2 = factor(tmp6$col_lab2, levels = niveau)
  return(tmp6)
})


pl= list()
for(y in seq_along(dt.map)){
  # y=2

  pl[[y]] = tts+
    ggnewscale::new_scale_fill() +
    geom_sf(data=dt.map[[y]],aes(fill=col_lab2), col = "NA")+
    scale_fill_manual("", values = c("white", cols2), na.value = "white" ) +
    guides(fill = guide_legend(reverse = T,
                               order=2,
                               title=nom,
                               nrow = 4,
                               title.hjust = 0,
                               title.vust = 0.5, size= 15))+
    geom_sf(data=xx, fill=NA, col="black")+
    theme_map()+
    theme(legend.key = element_rect(color = "black"),
          legend.position = "top",
          legend.spacing.x = unit(0.5, "cm"),
          legend.spacing.y = unit(1.5, "cm"),
          legend.title.align = 0.9,
          # legend.key.size = unit(0.5,"cm"),
          # legend.justification = "left",
          legend.direction = "horizontal",
          legend.title = element_text(size=16,
                                      margin = margin(r = 1, unit = "cm") ),
          legend.text = element_text(size=16,
                                     margin = margin(r = 0.7, unit = "cm") ),
          plot.margin = unit(c(0, 0, 0,0), "cm")# top, right, bottom, left
    )

}


ggarrange(plotlist = pl[6:1],common.legend = T,
          labels="AUTO",
          legend = "top")

jpeg(filename =file.path(workdir, "density_map_scenario.jpg"),
     width = 14,height = 7, units = "in", res=600,pointsize = 20)
ggarrange(plotlist = pl[6:1],common.legend = T,
          legend = "top")
dev.off()


tiff(filename = file.path(workdir, "density_map_scenario.tiff"),
     width = 12,height = 7, units = "in", res=600,pointsize = 12)
ggarrange(plotlist = pl[6:1],common.legend = T,
          legend = "top")
dev.off()
