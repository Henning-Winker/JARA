

# MAKE JABBA Multiplots
library("png")
# Set Working directory file where to store the results
File = "C:/Work/Research/GitHub/JARA"

# Read assessment species information file
sp.assess = read.csv(paste0(File,"/jara.assessments.csv"))
print(data.frame(Assessment.list=sp.assess[,1]))

#-----------------------------------------------------
# Examples of joining JARA outputs into multiplots
#-----------------------------------------------------

# All plot will be save in the respective species assessment folders

#----------------------------------------------------------------        
# Produce IUCN Supplement Information "style" Plot
#----------------------------------------------------------------
DIMs=c(5,4) # Plot dimension
sp=8  # Select species  
Plot = c("IUCN_SoupfinShark")
DIMs=c(2*5,2*4)/1.5
# Choos plot types
plots = c("AbundanceData_","Fits_","Residuals_","PopTrend_","AnnualRate_","IUCNplot_")[c(2,5,4,6)]
j=1
runs  = 1 # run name
# layout the plots into a matrix  columns, by row
par(mar=rep(0,4),omi= c(0, 0, 0, 0),cex=1.3) # no margins
layout(matrix(1:4, ncol=2, byrow=TRUE))
  
for(j in 1:length(plots)){
  run = runs[1]
  spsel= sp.assess[sp,]
  assessment = spsel$assessment # assessment name
  
  # plot path
  get_plot = paste0(File,"/",assessment,"/output",run,"/",plots[j],assessment,".png")
  # load image
  img <- readPNG(paste0(get_plot))
  
  # do the plotting
  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
  rasterImage(img,0,0,1,1)
  legend("topleft",paste0("(",letters[j],")"),bty="n",cex=1.2,x.intersp = -0.55,y.intersp=-0.2)
}
  
# Write plot
dev.print(png,paste0(File,"/",assessment,"/",Plot,"_JARA.png"), width = DIMs[1], height = DIMs[2], res = 300, units = "in")


#----------------------------------------------------------------        
# Produce  1 x 2 figure for Afr. Penguing trends
#----------------------------------------------------------------

    sp=c(1) # Select species  
    Plot = c("PenguinTrends")
    DIMs=c(2*5,1*4)
    # Choos plot types
    plots = c("AbundanceData_","Fits_","Residuals_","PopTrend_","AnnualRate_","IUCNplot_")[c(2,4)]
    j=1 
    runs  = 1 # run name
    # layout the plots into a matrix  columns, by row
    par(mar=rep(0,4),omi= c(0, 0, 0, 0),cex=1.3) # no margins
    layout(matrix(1:2, ncol=2, byrow=TRUE))
    
    for(j in 1:length(plots)){
        run = runs[1]
        spsel= sp.assess[sp,]
        assessment = spsel$assessment # assessment name
        
        # plot path
        get_plot = paste0(File,"/",assessment,"/output",run,"/",plots[j],assessment,".png")
        # load image
        img <- readPNG(paste0(get_plot))
        
        # do the plotting
        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
        rasterImage(img,0,0,1,1)
        #legend('topleft',paste0("a)"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
        #mtext("St. Joshep",side=3, outer=T, at=0.55,line=-1.2,cex=0.8)
        legend("topleft",paste0("(",letters[j],")"),bty="n",cex=1.2,x.intersp = -0.55,y.intersp=-0.2)
        }
    
    # Write plot
    dev.print(png,paste0(File,"/",assessment,"/",Plot,"_JARA.png"), width = DIMs[1], height = DIMs[2], res = 300, units = "in")

    
    #----------------------------------------------------------------        
    # Produce  1 x 2 figure for smoothhound shark projections
    #----------------------------------------------------------------
    
    p=c(7) # Select species  
    Plot = c("SmoothhoundPrj")
    DIMs=c(2*5,1*4)
    # Choos plot types
    plots = c("AbundanceData_","Fits_","Residuals_","PopTrend_","AnnualRate_","IUCNplot_")[c(2,4)]
    j=1 
    runs  = 1 # run name
    # layout the plots into a matrix  columns, by row
    par(mar=rep(0,4),omi= c(0, 0, 0, 0),cex=1.3) # no margins
    layout(matrix(1:2, ncol=2, byrow=TRUE))
    
    for(j in 1:length(plots)){
      run = runs[1]
      spsel= sp.assess[sp,]
      assessment = spsel$assessment # assessment name
      
      # plot path
      get_plot = paste0(File,"/",assessment,"/output",run,"/",plots[j],assessment,".png")
      # load image
      img <- readPNG(paste0(get_plot))
      
      # do the plotting
      plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
      rasterImage(img,0,0,1,1)
      #legend('topleft',paste0("a)"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
      #mtext("St. Joshep",side=3, outer=T, at=0.55,line=-1.2,cex=0.8)
      legend("topleft",paste0("(",letters[j],")"),bty="n",cex=1.2,x.intersp = -0.55,y.intersp=-0.2)
    }
    
    # Write plot
    dev.print(png,paste0(File,"/",assessment,"/",Plot,"_JARA.png"), width = DIMs[1], height = DIMs[2], res = 300, units = "in")
    
    
#----------------------------------------------------------------        
# Produce  1 x 2 figure for Afr. Penguing IUCN Status 
#----------------------------------------------------------------
    
    sp=c(1) # Select species  
    Plot = c("PenguinStatus")
    DIMs=c(2*5,1*4)
    # Choos plot types
    plots = c("AbundanceData_","Fits_","Residuals_","PopTrend_","AnnualRate_","IUCNplot_")[c(5:6)]
    j=1 
    runs  = 1 # run name
    # layout the plots into a matrix  columns, by row
    par(mar=rep(0,4),omi= c(0, 0, 0, 0),cex=1.3) # no margins
    layout(matrix(1:2, ncol=2, byrow=TRUE))
    
    for(j in 1:length(plots)){
      run = runs[1]
      spsel= sp.assess[sp,]
      assessment = spsel$assessment # assessment name
      
      # plot path
      get_plot = paste0(File,"/",assessment,"/output",run,"/",plots[j],assessment,".png")
      # load image
      img <- readPNG(paste0(get_plot))
      
      # do the plotting
      plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
      rasterImage(img,0,0,1,1)
      #legend('topleft',paste0("a)"),bty="n",x.intersp = -0.5,y.intersp = -0.2,cex=1.1)
      #mtext("St. Joshep",side=3, outer=T, at=0.55,line=-1.2,cex=0.8)
      legend("topleft",paste0("(",letters[j],")"),bty="n",cex=1.2,x.intersp = -0.55,y.intersp=-0.2)
    }
    
    # Write plot
    dev.print(png,paste0(File,"/",assessment,"/",Plot,"_JARA.png"), width = DIMs[1], height = DIMs[2], res = 300, units = "in")

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>    
# End of Code
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>    

    
    
    
    