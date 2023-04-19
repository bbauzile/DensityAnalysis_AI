#=============================================================================
# 			Density analysis paper of palmipeds farms in France
#
#  Mapping R0 in South west France
#  created by : Billy Bauzile (billy.bauzile@gmail.com)
#  originally created during PhD: 2019-2022
#
#
#=============================================================================



rm(list = ls())


library(sf)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(rgdal)
library(gstat)
library(raster)
# workdir
workdir = file.path("R0/")

#department
dep_fr=st_read("Data/Map_data/shapefile_national_mapFR.shp")#


sw = c("33", "40","64","65","32", "47", "24", "46", "82","81","24", "31",
       "11", "12", "34","48", "09", "16","17", "19", "15", "30")


com_fr= st_read(file.path("Data/Map_data/shapefile_southwFR_map.shp"))


dep_fr2 =  dep_fr |>
  filter(CODE_DEPT %in% sw)

com_fr2 = com_fr |>
  filter(CODE_DEPT %in% sw)


#
mp_dt2 = st_as_sf(dep_fr2) #just to change the name
mp_sw = st_as_sf(com_fr2) # just to change the variable name


# set epi vs non area

cnames <- c("id", "department", "insee", "type", "x", "y", "date_detection",
            "date_culling", "symptomatic", "preventively_culled", "delay_empty",
            "op", "is_pag", "in_zrp")

data=read.delim(file.path(workdir, "../Data/data.txt"),
                header = F,sep = " ", col.names = cnames)


mp_sw2 = merge(mp_sw, data |>
                 filter(date_detection>=0),
               by.x="INSEE_COM", by.y="insee", all=T)

mp_sw2v = st_as_sf(mp_sw2)

# Epizootic area

# The epidemic zone was created to account for the data implemented in [@Andronico2019].
# It containts farms in department that were in 100 km away from a detected case (galliformes or palmipeds farms).

geom= st_geometry(mp_sw2v[!is.na(mp_sw2v$id), ] )

x=st_union(st_buffer(geom, 100000), by_feature = F) #100km
xx= st_intersection(x,  mp_dt2, sparse = T) #epi area


mp_dt2$epi = 1

base_plot =ggplot()+
  geom_sf(data=mp_dt2,aes(fill=factor(epi) ))+
  geom_sf(data=xx, fill="grey95") +
  scale_fill_manual(values = c("#2B8CBE","grey95"),
                    name = "", #"Epidemic area",
                    labels = c("Outside\nepizootic area"),
                    guide = guide_legend(nrow=2, byrow = T))+
  theme(legend.text = element_text(size=16))
base_plot

#the R0
## data from the computation of R0
dt_r0 = list.files(path="Data", pattern = "niter",full.names = T)
nom = c("0.00282", "0.00345", "0.0043", "0.0071","0.01", "baseline")


data_r0= lapply(dt_r0, read_csv)
# data_r0[[1]]$sim |> max()


data_mod = list()

for(i in seq_along(data_r0)){
  data_mod[[i]] = data_r0[[i]] |>
    dplyr::filter(!is.na(R0), delta==11) |>
    left_join(data |>
                dplyr::select(id, insee, x,y),
              by=c("Farm"="id")) |>
    dplyr::group_by(scenario, delta, Farm) |>
    dplyr:: mutate(R_farm=quantile(R0, probs = 0.5, na.rm = T) ) |>
    distinct(Farm, .keep_all=T) |>
    group_by(insee) |>
    dplyr::mutate(R_com=quantile(R_farm, probs = 0.5) ) %>%
    distinct(insee, .keep_all=T) |>
    mutate(insee= as.character(insee),
           x= x*1000,
           y=y*1000) |>
    mutate(threshold =
             case_when(scenario=="baseline"~"Baseline",
                       scenario=='0.01'~"Reduction \n density to the 98th percentile",
                       scenario=="0.0071"~"Reduction \n density to the 95th percentile",
                       scenario=="0.0043"~"Reduction \n density to the 90th percentile",
                       scenario=="0.00345"~"Reduction \n density to the 85th percentile",
                       scenario=="0.00282"~"Reduction \n density to the 80th percentile"))

}


#sf
frc = st_transform(dep_fr2, CRS("+init=EPSG:2154"))

# Define the ouput grid - 50000 points in polygons extent
ptsreg = sp::spsample(as(xx, "Spatial"), 50000, type="regular")  |>
  st_as_sf()

f.1 = as.formula(R_farm~1)

# f.1  =list()
spdf = list()
# dat.fit = list()
Krig_ord = list()
Krig2 = list()
Krig_df = list()
for(d in seq_along(data_mod) ){
  # d=1
  #

  spdf[[d]] = SpatialPointsDataFrame(coords = data_mod[[d]][,c(7,8)],
                                    data= data_mod[[d]][,-c(7,8)])|>
    st_as_sf() |>
    st_set_crs('2154')


  st_crs(spdf[[d]])= st_crs(dep_fr2)
  st_crs(ptsreg)= st_crs(spdf[[d]])

  var.smpl <- variogram(f.1, spdf[[d]], cutoff=100000, width=500)

  # Compute the variogram model by passing the nugget, sill and range values
  # to fit.variogram() via the vgm() function.

  dat.fit <- fit.variogram(var.smpl, fit.ranges = T, fit.sills = T,
                           vgm(psill=0.4, model=c("Exp"), range=100000,nugget = 0.2 ))

  # The following plot allows us to assess the fit
  # plot(var.smpl, dat.fit, xlim=c(0,200000))
  Krig_ord[[d]] = krige(f.1, spdf[[d]], ptsreg, dat.fit,maxdist=15000)

  st_crs(Krig_ord[[d]]) = st_crs(dep_fr2)

  Krig2[[d]] = Krig_ord[[d]][!is.na(st_intersects(Krig_ord[[d]], dep_fr2)),]

}


p=list()
pal3 = rev(c('#fc8d59','#ffffbf','#91bfdb'))

r=list()
r.m = list()
map_prod = list()
for(dd in seq_along(Krig2)){
  #gridded only works for sp so there reconvert it
  Krig2[[dd]] = as(object = Krig2[[dd]], "Spatial")
  gridded(Krig2[[dd]])=T

  r[[dd]] <- raster(Krig2[[dd]])
  r.m[[dd]] <- mask(r[[dd]], dep_fr2)
  map_prod[[dd]] <- rasterToPoints(r.m[[dd]] ) |> as.data.frame() |>
    mutate(threshold =  unique(data_mod[[dd]]$threshold))
}


map_prod2=  do.call(rbind, map_prod)


colnames(map_prod2) = c("X","Y", "R0", "threshold")

vuni = unique(map_prod2$threshold)
vuni = vuni[6:1]
map_prod2$threshold = factor(map_prod2$threshold, levels = vuni)

## #########

map_prod2$R02 = cut(map_prod2$R0, breaks=c(0, 0.5, 1, 1.5, 2, 2.5) )

base_plot =ggplot()+
    geom_sf(data=mp_dt2 |> filter(!CODE_DEPT %in% c("07","66")),aes(fill=factor(epi) ))+
    geom_sf(data=xx, fill="#b2e2e2") +
    scale_fill_manual(values = c("#2B8CBE","#b2e2e2"),
                      name = "", #"Epidemic area",
                      labels = c("Outside\nepidemic area"),
                      guide = guide_legend(nrow=2, byrow = T))+
    theme(legend.text = element_text(size=14,
                                     margin = margin(r = 2, unit = "cm") ))+
    cowplot::theme_map()



#
# pal1 = c("#3182bd", "#deebf7","#feb24c", "#f03b20")
# pal3 = c("",'#91bfdb','#ffffbf','#fc8d59')

pal4 = c("#b2e2e2","#ffffb2", "#fed98e", "#fe9929")


lab = c("Reduction \n density to the 80th percentile"="F",
        "Reduction \n density to the 85th percentile"= "E",
        "Reduction \n density to the 90th percentile"="D",
        "Reduction \n density to the 95th percentile" ="C",
        "Reduction \n density to the 98th percentile" ="B",
        "Baseline"= "A")

# Figure used for density the article

R0 = expression(R[0])
## figure article ##
map_prodV2 = split(map_prod2, f = map_prod2$threshold)

carte = list()
for(i in seq_along(map_prodV2 ) ){
  carte[[i]] = base_plot+
    ggnewscale::new_scale_fill() +
    geom_raster(map_prodV2[[i]],
                mapping=aes(X,Y,
                            fill=R02 ))+
    facet_wrap("threshold~.", labeller = labeller(threshold=lab))+
    scale_fill_manual(R0,values = pal4, na.value="#fef0d9")+
    guides(fill=guide_legend(nrow = 2, reverse = T))+
    geom_sf(data=xx, fill=NA, col="black")+
    theme_map()+
    theme(legend.key = element_rect(color = "black"),
          legend.position = "top",
          legend.spacing.x = unit(0.5, "cm"),
          legend.spacing.y = unit(1.5, "cm"),
          legend.title.align = 0.9,
          legend.title = element_text(size=16,
                                      margin = margin(r = 0.5,
                                                      unit = "cm") ),
          legend.text = element_text(size=16,
                                     margin = margin(r = 2,
                                                     unit = "cm") ),
          legend.direction = "horizontal",
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0, 0,0), "cm")# top, right, bottom, left
    )

}


ggarrange(plotlist = carte[1:6], labels = "AUTO",
          font.label = list(size= 16), common.legend = T)

ggsave(filename = file.path("R0/figure_R0.png"), bg = "white")
