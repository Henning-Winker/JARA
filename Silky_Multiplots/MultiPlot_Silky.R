

# MAKE JABBA Multiplots
library("png")
# Set Working directory file where to store the results
File= "C:/Work/Research/Github/JARA"
# Set Assessment

setwd(File)  
multi.folder = "Silky_Multiplots"
dir.create(paste0(multi.folder),showWarnings = FALSE)
getwd()
assessment = c("Silky_Atl","Silky_Pac","Silky_Global")
  for(j in 1:3){
  Plot = assessment[j]
  DIMs=c(6,5)
  
  # setup plot
  par(mar=rep(0,4),omi= c(0, 0, 0, 0)) # no margins
  
  # layout the plots into a matrix w/ 12 columns, by row
  layout(matrix(1:4, ncol=2, byrow=TRUE))
  get_plot = paste0(Plot,"/Output/Fits_",Plot)
  
    # example image
    img <- readPNG(paste0(get_plot,".png"))
    
    # do the plotting
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
    legend('topleft',paste0("a)"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
    
    # Fit 
    get_plot = paste0(Plot,"/Output/lograte_",Plot)
    
    # example image
    img <- readPNG(paste0(get_plot,".png"))
    
    # do the plotting
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
    legend('topleft',paste0("b)"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
    k=1
    # Projection 
    scenario = paste0(Plot)  
    get_plot = paste0(scenario,"/Output/Poptrend_",scenario)
    
    # example image
    img <- readPNG(paste0(get_plot,".png"))
    
    # do the plotting
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
    legend('topleft',paste0(letters[k*2+1],")"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
    # ICUN
    get_plot = paste0(scenario,"/Output/IUCNplot_",scenario)
    
    # example image
    img <- readPNG(paste0(get_plot,".png"))
    
    # do the plotting
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
    legend('topleft',paste0(letters[(k+1)*2],")"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
    
  
  dev.print(png,paste0(multi.folder,"/",Plot,"_JARA.png"), width = DIMs[1], height = DIMs[2], res = 200, units = "in")
}







  