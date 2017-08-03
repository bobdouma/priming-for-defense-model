# script to make simulations on the benefits of priming

# remove objects from workspace
rm(list=ls(all=TRUE))

library(reshape2)
library(gridExtra)
library(grid)



# load model
source("D:/BobDouma/Veni/Game_theory_plantcommunication/R-scripts/20160301_Flower_Plant_competition_model_defense_v5_clean.R")
setwd("D:/BobDouma/Veni/Game_theory_plantcommunication/simulations/sim_8_linear_delay6_lognormal_fixC_def_Ft")
#source("M:/My Documents/Plantcommunication/20160301_Flower_Plant_competition_model_defense_v5_clean.R")
#setwd("M:/My Documents/Plantcommunication/simulations/sim_8_linear_delay6_lognormal_fixC_def_Ft")


##### Analysis 1

# characteristics:
# fixed flowering time
# natural pop dynamics for feeding rate
# C_def max = 8
 
# set parameters
fPrim = c(0)
C_def = 1 # max is one
D = c(1,2,6,8)
Ths = c(10,20)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
#Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Co_r = c(6*1e-08,22*1e-08)
Thd = 6
Tpd = 100
alpha = c(1,1.25)#c(1,1.1,1.25)
feeding = c("lognorm")
responsepat = c("Fa","Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("20160412_Input_defense_evo_Priming_onoff.csv")
pars = as.data.frame(matrix(NA,nrow=nrow(a),ncol=length(comb[1,])))

names(pars) = names(comb)
pars[] = comb[1,]

pars$D = a$D
pars$fPrim_f = a$fPrim
pars$fPrim_r = a$fPrim
pars$C_prim_f = a$C_prim
pars$C_prim_r = a$C_prim
pars$Ft = a$Ft
pars$mod = "Prim.f"
pars$Th_s = a$Th_s
pars$Th_e = a$Th_e
pars$Co_r = a$Co_r
pars$Tp_s = a$Tp_s
pars$Tp_ss = a$Tp_ss
pars$Tp_dd = a$Tp_dd
pars$C_def = a$C_def
pars$C_def_f = a$C_def
pars$C_def_r = a$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = 6
pars$Tco_d = 12
pars$shape = "delay"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat

# if Co_r is 0, Def.f and Def.r FALSE
temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1.1e-07,]
pars$ref_Co_r = pars$Co_r
temp$ref_Co_r = temp$Co_r
temp$Co_r = 0
temp[,c("Def.f","Def.r","H.f","H.r")] = FALSE
temp$no = seq(from=(max(pars$no)+1),to=(nrow(temp)+max(pars$no)))
pars = rbind(pars,temp)
# optimal flowering time without herbivory
Opt.flower = data.frame(D=c(1,2,6,8),Ft=c(90.58333333,84.875,75.83333333,73.45833333))
pars$Ft = Opt.flower[match(pars$D,Opt.flower$D),2]
pars$Ft = 80
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
# add outcomes of simulations to parameters
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)

# store results
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft.csv",sep=""))

#out.scenario = read.csv("20170726_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft.csv")

# relative benefit of priming compared to non-priming
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100

# priming response faster
dat= data.frame(y=0)
out.scenario.fa = out.scenario[out.scenario$responsepat == "Fa",]
pdf("Priming_beneficial_attack_linear_8_delay6_Cor_12d_scenarios_lognormal_fixC_def_Ft_fa.pdf")
ggplot() +
  geom_line(aes(x=C_prim_f,y=rel,colour=factor(D)),data=out.scenario.fa) +
  facet_grid(Co_r+alpha~Th_s+Tp_ss,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9)) +
  ylab("Relative change (%)")
dev.off()

# priming response earlier
out.scenario.ea = out.scenario[out.scenario$responsepat == "Ea",]
pdf("Priming_beneficial_attack_linear_8_delay6_Cor_12d_scenarios_lognormal_fixC_def_Ft_ea.pdf")
ggplot() +
  geom_line(aes(x=C_prim_f,y=rel,colour=factor(D)),data=out.scenario.ea) +
  facet_grid(Co_r+alpha~Th_s+Tp_ss,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9)) +
  ylab("Relative change (%)")
dev.off()

