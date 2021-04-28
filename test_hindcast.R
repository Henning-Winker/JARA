#-----------------------------------------
# Test JARA hindcast options
#-----------------------------------------

# UPDATE JARA BEFORE TESTING
# devtools::install_github("henning-winker/JARA")

library(JARA)
data("jaradata")

output.dir = "C:/Work/Research/GitHub/JARA_Examples/Hindcast"
# Run hindcasts for yellowfin skate
d = jrdat$Yellowspot_skate
yps = jarainput= build_jara(I=d$I,se=d$SE,assessment="YellowspottedSkate",
                 GL=12)
fit1 = fit_jara(yps,output.dir = output.dir)
jrplot_fits(fit1)
jrplot_trjfit(fit1)
jrplot_poptrj(fit1)
jrplot_state(fit1) # new plot
# save 
jrplot_state(fit1,output.dir = output.dir,as.png = T)

# use inbuilt hindcast function (refit retrospectively)
hc1 = jara_hindcast(yps,plotall = T,output.dir = output.dir,peels=0:5,speedup = F)
# try out plot
jrplot_retrobias(hc1) # Retrospective bias
jrplot_retroiucn(hc1) # Reptropective threat perception
# save
jrplot_retrobias(hc1,as.png = T,output.dir = output.dir) # Retrospective bias
jrplot_retroiucn(hc1,as.png = T,output.dir = output.dir) # Reptropective threat perception

# Run hindcasts for yellowfin skate
d = jrdat$SmoothhoundShark
shs = build_jara(I=d$I,se=d$SE,assessment="Smoothhound",GL=13.7)
fit2 = fit_jara(shs,output.dir = output.dir)
jrplot_fits(fit2)
jrplot_trjfit(fit2)
jrplot_poptrj(fit2)
jrplot_state(fit2) # new plot
# save 
jrplot_state(fit2,output.dir = output.dir,as.png = T)


# Run on fast mode (Half GT)
hc2 = jara_hindcast(shs,plotall = T,output.dir = output.dir,peels=0:5,speedup = T)
# try out plot
jrplot_retrobias(hc2) # Retrospective bias
jrplot_retroiucn(hc2) # Reptropective threat perception
# save
jrplot_retrobias(hc2,as.png = T,output.dir = output.dir) # Retrospective bias
jrplot_retroiucn(hc2,as.png = T,output.dir = output.dir) # Reptropective threat perception


