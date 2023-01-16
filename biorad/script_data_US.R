
library(bioRad)
# library(bioRad, lib.loc="C:/Users/rnussba1/Documents/GitHub/bioRad")

library(ggplot2)
library(tidyr)
library('plot.matrix')
library(profvis)

# folder <- 'C:/Users/rnussba1/Box/BMM-US'
folder <- '/Users/raphael/Library/CloudStorage/Box-Box/BMM-US/'
setwd(folder)

files <- list.files('../nexrad_vpts/')

radar.list <- unique(unlist(lapply(files, function(x){substring(x, 1, 4)})))
year.list <- unique(unlist(lapply(files, function(x){substring(x, 5, 8)})))
# year.list <- c('2019')#

ex <- matrix(nrow =length(radar.list), ncol = length(year.list), dimnames = list(radar.list, year.list))
lat <- c()
lon <- c()
height <- c()
vpts.prof <- matrix(nrow =length(radar.list), ncol = 50 )

time_full <- seq(as.POSIXct('1995-01-01', tz="UTC"),as.POSIXct('2021-12-31 23:45', tz="UTC"), by=15*60)
vid <- matrix(nrow =length(radar.list), ncol = length(time_full), dimnames = list(radar.list, format(time_full)))
viu <- matrix(nrow =length(radar.list), ncol = length(time_full), dimnames = list(radar.list, format(time_full)))
viv <- matrix(nrow =length(radar.list), ncol = length(time_full), dimnames = list(radar.list, format(time_full)))
viwu <- matrix(nrow =length(radar.list), ncol = length(time_full), dimnames = list(radar.list, format(time_full)))
viwv <- matrix(nrow =length(radar.list), ncol = length(time_full), dimnames = list(radar.list, format(time_full)))

if (FALSE){
  vid <- t(read.csv('data/vid.csv'))
  vid <- cbind(vid, matrix(nrow =length(radar.list), ncol = length(time_full)-dim(vid)[2]))
  viv <- t(read.csv('data/viv.csv'))
  viv <- cbind(viv, matrix(nrow =length(radar.list), ncol = length(time_full)-dim(viv)[2]))
  viu <- t(read.csv('data/viu.csv'))
  viu <- cbind(viu, matrix(nrow =length(radar.list), ncol = length(time_full)-dim(viu)[2]))
  viwv <- t(read.csv('data/viwv.csv'))
  viwv <- cbind(viwv, matrix(nrow =length(radar.list), ncol = length(time_full)-dim(viwv)[2]))
  viwu <- t(read.csv('data/viwu.csv'))
  viwu <- cbind(viwu, matrix(nrow =length(radar.list), ncol = length(time_full)-dim(viwu)[2]))
}


#i_y=which(year.list==2018)
# i_r = which(radar.list=="KBBX")

for (i_y in 1:length(year.list)){
  
  index_time <- format(time_full, format="%Y")==year.list[i_y]
  
  for (i_r in 1:length(radar.list)){
  
    filename <- paste0('../nexrad_vpts/', radar.list[i_r], year.list[i_y] ,'.rds')
    ex[i_r,i_y] <- file.exists(filename)
    
    
    if (ex[i_r,i_y]){
      print(paste0(i_r,' ',filename))
      vpts <- readRDS(filename)
      
      # plot(vpts[25000+(1:700)])
      
      #lat[i_r]=vpts$attributes$where$lat
      #lon[i_r]=vpts$attributes$where$lon
      #height[i_r]=vpts$attributes$where$height

      sd_vvp_threshold(vpts) <- 0

      vpts.reg <- regularize_vpts(vpts,
                                  date_min=as.POSIXct(paste0(year.list[i_y],'-01-01'), tz="UTC"),
                                  date_max=as.POSIXct(paste0(year.list[i_y],'-12-31 23:45'), tz="UTC"),
                                  interval = 15,
                                  fill = 7.5,
                                  units="mins")
      
      vpts.reg.f = filter_dbzh(vpts.reg, threshold=500, height=2000, agl_max=Inf, drop=F, quantity="dens")
      
      
      # vpts.prof[i_r,] = rowSums(vpts.reg.f$data$dens,na.rm=T)
      
      vpits <-integrate_profile(vpts.reg.f, alt_min="antenna")
      
      # u=u+4*24*8
      # 
      # m=100
      # par(mfrow=c(4,1), mar = c(2, 2, 0, 0))
      # span = u+(-m:m);
      # plot(vpts.reg[span])
      # plot(vpts.reg[span],quantity = 'DBZH')
      # plot(vpts.reg.f[span])
      # plot(vpits$datetime[span] ,vpits$vid[span],ylab="dens integrated",xaxs="i")
      # 
      # #
      vid[i_r,index_time] = vpits$vid
      viu[i_r,index_time] = vpits$u
      viv[i_r,index_time] = vpits$v
      viwu[i_r,index_time] = vpits$u_wind
      viwv[i_r,index_time] = vpits$v_wind
      #viwv[i_r,index_time] = vpits$sd_vvp
      
      # tmp <- t(vpts$data$eta)
      # rownames(tmp) <- vpts$datetime
      # colnames(tmp) <- vpts$height
      # tmp <- tmp[,rowSums(t(is.na(tmp)))!=dim(tmp)[1]]
      # write.csv(tmp,paste0('data/vpts-csv/',radar.list[i_r], year.list[i_y],'-eta.csv'))
      # 
      # tmp <- t(vpts$data$sd_vvp)
      # rownames(tmp) <- vpts$datetime
      # colnames(tmp) <- vpts$height
      # tmp <- tmp[,rowSums(t(is.na(tmp)))!=dim(tmp)[1]]
      # write.csv(tmp,paste0('data/vpts-csv/',radar.list[i_r], year.list[i_y],'-sdvvp.csv'))
      
      # tmp <- t(vpts$data$airspeed)
      # rownames(tmp) <- vpts$datetime
      # colnames(tmp) <- vpts$height
      # tmp <- tmp[,rowSums(t(is.na(tmp)))!=dim(tmp)[1]]
      # write.csv(tmp,paste0('data/vpts-csv/',radar.list[i_r], year.list[i_y],'-airspeed.csv'))
    }
  }
  
  # save(vid,viv,viu,viwu,viwv,file="viduv.Rdata")
  print(paste0('Saved ',i_y))
}

write.csv(t(vid),'data/vid.csv',row.names = FALSE,col.names = FALSE)
write.csv(t(viv),'data/viv.csv',row.names = FALSE,col.names = FALSE)
write.csv(t(viu),'data/viu.csv',row.names = FALSE,col.names = FALSE)
write.csv(t(viwv),'data/viwv.csv',row.names = FALSE,col.names = FALSE)
write.csv(t(viwu),'data/viwu.csv',row.names = FALSE,col.names = FALSE)
plot(ex)




profvis({
  
  vpts <- readRDS(filename)
  vpts.reg <- regularize_vpts(vpts, interval = 15, units="mins")
  
  
  
  vpits <-integrate_profile(vpts.reg)

})

plot(vpts$height,colSums(vpts.prof,na.rm = T))

vpts.exp <- vpts.reg
attr(vpts.exp, "class") <- NULL
vpts.exp$timesteps <-c()
x <- toJSON(vpts.exp)
write_json(x,'KABR2020.json')


filter_dbzh <- function(vpts, threshold=7,height=1000, agl_max=Inf, drop=F, quantity="DBZH"){
  height_index_max <- ((vpts$attributes$where$height + agl_max) %/% vpts$attributes$where$interval)
  height_index_max <- min(dim(vpts)[2],height_index_max)
  height_range <- colSums(vpts$data[[quantity]][1:height_index_max,]>threshold,na.rm=T)*vpts$attributes$where$interval
  index <- which(height_range > height)
  if(length(index)==0) return(vpts)
  # if remove, drop the profiles
  if(drop) return(vpts[-index])
  # otherwise set the density field to NA, but keep the profile
  vpts$data$DBZH[,index] <- NA
  vpts$data$dens[,index] <- NA
  vpts
}

