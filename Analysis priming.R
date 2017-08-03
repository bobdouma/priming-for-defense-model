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

# analyse difference of faster response to earlier response
out.scenario.subset = out.scenario[out.scenario$Th_s ==10 & out.scenario$C_prim_f ==0.5 & out.scenario$D == 8,]
out.scenario.subset = out.scenario[]
out.responsepat = out.scenario.subset[,c("rel","responsepat")] 

out.responsepat =out.responsepat[out.responsepat$rel < 200,]
plot(out.responsepat[out.responsepat$responsepat =="Ea","rel"]~ out.responsepat[out.responsepat$responsepat =="Fa","rel"])
summary(lm(out.responsepat[out.responsepat$responsepat =="Ea","rel"]~ out.responsepat[out.responsepat$responsepat =="Fa","rel"]))
abline(lm(out.responsepat[out.responsepat$responsepat =="Ea","rel"]~ out.responsepat[out.responsepat$responsepat =="Fa","rel"]))


############ apply error management theory to results to plant density
##############
# set parameters
fPrim = c(0)
C_def = 8
D = c(1,2,3,4,5,6,7,8)
Ths = c(10,15,20)
The = c(180)
Tps = c(-8)
C_prim =  c(0.05,0.1,0.25)
Co_r = c(6*1e-08,22*1e-08)
Thd = 6
Tpd = 8
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
temp = pars[pars$fPrim_f ==0 & pars$Co_r ==2.2e-07,]
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
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_EMT.csv",sep=""))


############ apply error management theory to results to plant age
##############
# set parameters
fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(5,6,7,8,910,11,12,13,14,15,16,17,18,19,20)
The = c(180)
Tps = c(-4,-8)
C_prim =  c(0.05,0.1,0.25)
Co_r = c(6*1e-08,22*1e-08)
Thd = 6
Tpd = 8
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
temp = pars[pars$fPrim_f ==0 & pars$Co_r ==2.2e-07,]
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
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_EMT.csv",sep=""))


########################################################################################################
### Simulations 2
### make herbivory proportional to LAI

# set parameters
fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(3,10,20,30)
The = c(180)
Tps = c(-2,-8,-12)
C_prim =  c(0,0.1,0.25,0.5,0.75,1)
Co_r = c(5e-08,1e-07,7.5e-07,5e-06)
Thd = 6
Tpd = 100
alpha = c(1,1.1,1.25)
feeding = c("prop")
responsepat = c("Fa","Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1e-07,]
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
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_propLAI_fixC_def_Ft.csv",sep=""))

out.scenario = read.csv("20160727_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_propLAI_fixC_def_Ft.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s

out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100

# priming response faster
dat= data.frame(y=0)
out.scenario.fa = out.scenario[out.scenario$responsepat == "Fa",]
pdf("Priming_beneficial_attack_linear_8_delay6_Cor_12d_scenarios_propLAI_fixC_def_Ft_fa.pdf")
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
pdf("Priming_beneficial_attack_linear_8_delay6_Cor_12d_scenarios_propLAI_fixC_def_Ft_ea.pdf")
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

# benefits greater when Th_s later?
bef = out.scenario.ea[out.scenario.ea$Co_r==5.0e-06 & out.scenario.ea$Th_s %in% c(20,30) & out.scenario.ea$Tp_ss == -2 & out.scenario.ea$D ==8 & out.scenario.ea$C_prim_f == 0.5,]
bef = out.scenario.ea[out.scenario.ea$Co_r==0 & out.scenario.ea$Th_s %in% c(20,30) & out.scenario.ea$Tp_ss == -2 & out.scenario.ea$D ==8 & out.scenario.ea$C_prim_f == 0.5,]

# relative benefits and relative costs are lower later in the season...
# benefits and costs of priming are higher in high density. Is this due to initial higher LAI?

pars.sub = pars[pars$Co_r %in% c(0,7.5e-07) & pars$Th_s==20 & pars$Tp_s ==8 & pars$C_prim_f ==0.5 & pars$responsepat == "Ea",]
pars.sub$S_m = pars.sub$S_m/pars.sub$D

results1 = find_endpoints_defense(pars.sub,Ft=F)
results1 = data.frame(pars.sub,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_propLAI_fixC_def_adj_Sm_Ft.csv",sep=""))

# no this is not due to higher initial LAI, but due to the fact that a abolsute reduction in dLAI (which is constant over densities),
# leads to a relative higher reduction in Pnet high density than in low density 


########################################################################################################
### Simulations 3
### When is priming more beneficial; given the same seasonal herbivory pressure; immediate herbivory, or a build up of herbivores?


Ts = 180
fPrim = c(1)
C_def = 8
D = c(1,2,6,8)
Ths = c(20)
The = c(100,180)
Tps = c(-2,-8,-12)
C_prim =  c(0,0.1,0.25,0.5,0.75,1)
Co_r = c(5e-08,1e-08)
Thd = 6
Tpd = 100
alpha = c(1)
Tco_d = c(50,25,12,6,0) 
feeding = c("constant")
responsepat = c("Fa","Ea")
# make parameter matrix
a = expand.grid(Ts=Ts,C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,Tco_d=Tco_d,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
#a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Co_r.or = a$Co_r
pars$Co_r = (a$Co_r*Ts)/(0.5*((a$Tco_d+a$Th_s)-Th_s.i)+(Th_e - (a$Tco_d+a$Th_s)))
pars$Tp_s = a$Tp_s
pars$Tp_ss = a$Tp_ss
pars$Tp_dd = a$Tp_dd
pars$C_def = a$C_def
pars$C_def_f = a$C_def
pars$C_def_r = a$C_def
pars$Ts = a$Ts
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = 6
pars$Tco_d = a$Tco_d
pars$shape = "delay"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat

# if Co_r is 0, Def.f and Def.r FALSE
# temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1.1e-07,]
# pars$ref_Co_r = pars$Co_r
# temp$ref_Co_r = temp$Co_r
# temp$Co_r = 0
# temp[,c("Def.f","Def.r","H.f","H.r")] = FALSE
# temp$no = seq(from=(max(pars$no)+1),to=(nrow(temp)+max(pars$no)))
# pars = rbind(pars,temp)
# # optimal flowering time without herbivory
# Opt.flower = data.frame(D=c(1,2,6,8),Ft=c(90.58333333,84.875,75.83333333,73.45833333))
# pars$Ft = Opt.flower[match(pars$D,Opt.flower$D),2]
pars$Ft = 80
# do simulations
# subset
#pars.sub = pars[pars$Th_s ==20 & pars$Co_r.or == 1e-08 & pars$D==8 & pars$Tp_s == 12 & pars$responsepat == "Ea" & pars$Tco_d %in% c(0,6,50),]
#results1 = find_endpoints_defense(pars.sub,Ft=F)
#results1 = data.frame(pars.sub,results1)
#results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results

write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_varying_Tco_d.csv",sep=""))

out.scenario = read.csv("20160728_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_varying_Tco_d.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s

out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100

# priming response earlier
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100
out.scenario.ea = out.scenario[out.scenario$responsepat == "Ea",]
ggplot() +
  geom_line(aes(x=C_prim_f,y=-seeds_tfm,colour=factor(D)),data=out.scenario.ea) +
  facet_grid(Tco_d+Co_r.or ~Th_s+Th_e+Tp_ss,scales="free") +
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


###########################
### Simulations 4
###### response of priming faster, earlier and stronger?

# set parameters
fPrim = c(0)
C_def = c(4,8)
D = c(1,8)
Ths = c(10)
The = c(180)
Tps = c(-7)
C_prim =  c(0,0.125,0.25,0.375,0.5,0.675,0.75,0.875,1)
Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Thd = 6
Tpd = 100
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Fa","Ea","St")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
#a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
pars = as.data.frame(matrix(NA,nrow=nrow(a),ncol=length(comb[1,])))

names(pars) = names(comb)

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
# temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1.1e-07,]
# pars$ref_Co_r = pars$Co_r
# temp$ref_Co_r = temp$Co_r
# temp$Co_r = 0
# temp[,c("Def.f","Def.r","H.f","H.r")] = FALSE
# temp$no = seq(from=(max(pars$no)+1),to=(nrow(temp)+max(pars$no)))
# pars = rbind(pars,temp)
# optimal flowering time without herbivory
Opt.flower = data.frame(D=c(1,2,6,8),Ft=c(90.58333333,84.875,75.83333333,73.45833333))
pars$Ft = Opt.flower[match(pars$D,Opt.flower$D),2]
pars$Ft = 80

# adjust simulations for "Faster"
pars[pars$responsepat == "St","fPrim_f"] = pars[pars$responsepat == "St","C_prim_f"]

# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_ResponsePattern_v2.csv",sep=""))


out.scenario = read.csv("20160728_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_ResponsePattern_v2.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s

out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100

# priming response earlier
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100
out.scenario = out.scenario[out.scenario$D == 8,]

dat= data.frame(y=0)
ggplot() +
  geom_line(aes(x=C_prim_f,y=-seeds_tfm,colour=factor(responsepat)),data=out.scenario) +
  facet_grid(Co_r+alpha ~ C_def,scales="free") +
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


ggplot() +
  geom_line(aes(x=C_prim_f,y=seeds_tfm,colour=factor(responsepat)),data=out.scenario) +
  facet_grid(Co_r+alpha ~Th_s+Th_e+Tp_ss,scales="free") +
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


out.scenario = read.csv("20160802_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_ResponsePattern_v2.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100
out.scenario = out.scenario[out.scenario$D ==8,]
out.scenario = out.scenario[out.scenario$alpha ==1,]
out.scenario$Co_r = factor(out.scenario$Co_r)
pdf("Response_patterns.pdf")
ggplot() +
  geom_line(aes(x=C_prim_f,y=-seeds_tfm,colour=Co_r),data=out.scenario) +
  facet_grid(responsepat ~C_def,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), , strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9)) +
  ylab("difference in fitness between primed and non-primed focal plant")
dev.off()

#pars.subset = pars[pars$Co_r == 2.2e-07 & pars$responsepat %in% c("Ea","St","Fa") & pars$D==8 & pars$Th_s == 30 & pars$alpha ==1 & pars$C_prim_f %in% c(0,1) & pars$Tp_ss == -2,]
test = c(757,758,759,760,2485,2486,2487,2488,4213,4214,4215,4216)

pars.subset = pars[test,]
pars.subset$Th_d = 6
pars.subset$Tp_s =8
pars.subset$Th_s =10
pars.subset$feeding = "lognorm"
pars.subset$C_def = 7.5
pars.subset$C_def_f = pars.subset$C_def
pars.subset$C_def_r = pars.subset$C_def
pars.subset$Co_r = 2e-08
results2 = find_endpoints_defense(pars.subset,Ft=F)
results2

 

############################################
########### Simulations 5 #################
# costs and benefits of priming when different priming durations

# set parameters
fPrim = c(0)
C_def = 8
D = c(1,8)
Ths = c(20,30)
The = c(180)
Tps = seq(from=-19,to=-1,length.out=10)
C_prim =  c(0,0.1,0.25,0.5,0.75,1)
Co_r = c(11*1e-08,22*1e-08)
Thd = 6
Tpd = 100
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Fa","Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
#a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1.1e-07 & pars$responsepat =="Ea",]
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
# adjust simulations for "Faster"
pars[pars$responsepat == "St","fPrim_f"] = pars[pars$responsepat == "St","C_prim_f"]

# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying.csv",sep=""))

out.scenario = read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s

out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100
out.scenario = out.scenario[out.scenario$alpha ==1 & out.scenario$D ==8,]

dat= data.frame(y=0)

ggplot() +
  geom_line(aes(x=C_prim_f,y=rel,colour=factor(responsepat)),data=out.scenario) +
  facet_grid(Co_r+responsepat ~Th_s+Th_e+Tp_ss,scales="free") +
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

out.scenario.sub = out.scenario[out.scenario$Co_r %in% c(0,1.1e-07) & out.scenario$Th_s == 20,]
subset = out.scenario.sub[out.scenario.sub$responsepat == "Ea" & out.scenario.sub$Co_r ==0,]
subset$responsepat = "Fa"
out.scenario.sub = rbind(out.scenario.sub,subset)
ggplot()+
  geom_line(aes(x=Tp_ss, y =rel, colour=factor(Co_r)), data=out.scenario.sub)+
  facet_grid(responsepat ~ C_prim_f )+
  geom_hline(aes(yintercept=y),data=dat) 
  
#########################
###### simulations 6


########### Simulations 5 #################
# costs and benefits of priming when different priming durations

# set parameters
fPrim = c(0)
C_def = 8
D = c(8)
Ths = seq(from=1,to=16,1)
The = c(180)
Tps = c(3,10)
C_prim =  seq(0,1,length.out = 50)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = c(6*1e-08,22*1e-08)
Thd = 6
Tpd = seq(from=1,to=16,1)
alpha = c(1)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Th_s = a$Th_s + a$Tp_s
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
#a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s,(a$Tp_d+a$Tp_s)) - a$Tp_s
#a = unique(a[,-15])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
temp = pars[pars$Co_r ==2.2e-07 & pars$responsepat =="Ea",]
pars$ref_Co_r = pars$Co_r
temp$ref_Co_r = temp$Co_r
temp$Co_r = 0
temp$C_prim_f = 0
temp$C_prim_r = 0
temp[,c("Def.f","Def.r","H.f","H.r")] = FALSE
temp$no = seq(from=(max(pars$no)+1),to=(nrow(temp)+max(pars$no)))
temp$Th_s = 50
# temp$Ft = 80
#pars =temp
pars = rbind(pars,temp)
# optimal flowering time without herbivory
Opt.flower = data.frame(D=c(1,2,6,8),Ft=c(90.58333333,84.875,75.83333333,73.45833333))
pars$Ft = Opt.flower[match(pars$D,Opt.flower$D),2]
pars$Ft = 80
# adjust simulations for "Faster"
pars[pars$responsepat == "St","fPrim_f"] = pars[pars$responsepat == "St","C_prim_f"]

pars[pars$Co_r ==0,"Tp_d"] = pars[pars$Co_r ==0,"Tp_dd"]
plot(pars$Tp_d ~jitter(pars$Th_s),col=ifelse(pars$Co_r>0,1,0)+1)
# do simulations
#results1 = find_endpoints_defense(pars,Ft=F)
results1 = find_endpoints_defense(temp,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying_optimal_Tpd_Cor0_Cprim0.csv",sep=""))


### analysis

out.scenario = read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying_optimal_Tpd.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100


out.scenario.mean = melt(tapply(out.scenario$seeds_mu,list("D"=out.scenario$D,
                                                       "Tp_d"=out.scenario$Tp_dd,
                                                       "Co_r"=out.scenario$Co_r,
                                                       "C_prim"=out.scenario$C_prim_f),mean))

out.scenario.mean.tf = melt(tapply(out.scenario$seeds_tf,list("D"=out.scenario$D,
                                                            "Tp_d"=out.scenario$Tp_dd,
                                                            "Co_r"=out.scenario$Co_r,
                                                            "C_prim"=out.scenario$C_prim_f),mean))

out.scenario.mean = data.frame(out.scenario.mean,out.scenario.mean.tf[,5])
names(out.scenario.mean)[5:6] = c("mu","tf")
out.scenario.mean$rel = (out.scenario.mean$mu-out.scenario.mean$tf)/(out.scenario.mean$tf)*100


out.scenario[out.scenario$D ==8 & out.scenario$C_prim_f == 0.5 & out.scenario$Co_r ==0,]

dat= data.frame(y=0)

ggplot() +
  geom_line(aes(x=Tp_d,y=rel,colour=factor(C_prim)),data=out.scenario.mean) +
  facet_grid(Co_r ~ D,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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


out.scenario.15 = read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying_optimal_Tpd.csv")
out.scenario.15$Tp_ss = out.scenario.15$Tp_s - out.scenario.15$Th_s
out.scenario.15$rel = -out.scenario.15$seeds_tfm/(out.scenario.15$seeds_tf)*100
out.scenario.15 = out.scenario.15[out.scenario.15$Th_s < 19,]
out.scenario.15 = out.scenario.15[out.scenario.15$Tp_dd < 17,]


out.scenario.mean.15 = melt(tapply(out.scenario.15$seeds_mu,list("D"=out.scenario.15$D,
                                                           "Tp_d"=out.scenario.15$Tp_dd,
                                                           "Co_r"=out.scenario.15$Co_r,
                                                           "C_prim"=out.scenario.15$C_prim_f),mean))

out.scenario.mean.tf.15 = melt(tapply(out.scenario.15$seeds_tf,list("D"=out.scenario.15$D,
                                                              "Tp_d"=out.scenario.15$Tp_dd,
                                                              "Co_r"=out.scenario.15$Co_r,
                                                              "C_prim"=out.scenario.15$C_prim_f),mean))

out.scenario.mean.15 = data.frame(out.scenario.mean.15,out.scenario.mean.tf.15[,5])
names(out.scenario.mean.15)[5:6] = c("mu","tf")
out.scenario.mean.15$rel = (out.scenario.mean.15$mu-out.scenario.mean.15$tf)/(out.scenario.mean.15$tf)*100
out.scenario.mean.15 = out.scenario.mean.15[out.scenario.mean.15$D==8,]

dat= data.frame(y=0)
ggplot() +
  geom_line(aes(x=Tp_d,y=value,colour=factor(C_prim)),data=out.scenario.mean.prob) +
  facet_grid(variable~.,scales="free") +
  #  geom_hline(aes(yintercept=y),data=dat) +
  geom_point(aes(x=Tp_d,y=value,colour=factor(C_prim)),data=max.vals)+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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



###### 15 days and x probability of no attack
prob = 0.3
out.scenario.mean.prob0 = dcast(out.scenario.mean.15[,c(1,2,3,4,7)], D+Tp_d+C_prim~Co_r)
out.scenario.mean.prob0$`6e-08 P(attack) = 0.7` = prob*out.scenario.mean.prob0$`0` + (1-prob)*out.scenario.mean.prob0$`6e-08`
out.scenario.mean.prob0$`0 P(attack) = 0` = 1*out.scenario.mean.prob0$`0` + 0*out.scenario.mean.prob0$`6e-08`
out.scenario.mean.prob0$`6e-08 P(attack) = 1` = 0*out.scenario.mean.prob0$`0` + (1)*out.scenario.mean.prob0$`6e-08`

out.scenario.mean.prob = melt(out.scenario.mean.prob0[,c(1,2,3,7,8,9)],id.var=1:3)
out.scenario.mean.prob$variable = relevel(out.scenario.mean.prob$variable,ref = 2)

out.scenario.mean.prob.max = melt(tapply(out.scenario.mean.prob$value,list(out.scenario.mean.prob$variable,out.scenario.mean.prob$D),which.max))
max.val = lapply(split(out.scenario.mean.prob,list(out.scenario.mean.prob$variable,out.scenario.mean.prob$D)),function(x){a = data.frame(max(x$value),x[which.max(x$value),2]);return(a)})
max.vals = data.frame(c("0 P(attack) = 0","6e-08 P(attack) = 0.7","6e-08 P(attack) = 1"),c(1,3,15),c(0,0.1519322,0.5767987))
names(max.vals) = c("variable","Tp_d","value")
max.vals$C_prim = 0.1
dat= data.frame(y=0)
ggplot() +
  geom_line(aes(x=Tp_d,y=value,colour=C_prim),data=out.scenario.mean.prob) +
  facet_grid(variable~.,scales="free") +
#  geom_hline(aes(yintercept=y),data=dat) +
  geom_point(aes(x=Tp_d,y=value),data=max.vals)+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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




#### 5 days

out.scenario.8 = read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying_optimal_Tpd.csv")
out.scenario.8$Tp_ss = out.scenario.8$Tp_s - out.scenario.8$Th_s
out.scenario.8$rel = -out.scenario.8$seeds_tfm/(out.scenario.8$seeds_tf)*100
out.scenario.8 = out.scenario.8[out.scenario.8$Th_s < 12,]
out.scenario.8 = out.scenario.8[out.scenario.8$Tp_dd < 8,]


out.scenario.mean.8 = melt(tapply(out.scenario.8$seeds_mu,list("D"=out.scenario.8$D,
                                                                 "Tp_d"=out.scenario.8$Tp_dd,
                                                                 "Co_r"=out.scenario.8$Co_r,
                                                                 "C_prim"=out.scenario.8$C_prim_f),mean))

out.scenario.mean.tf.8 = melt(tapply(out.scenario.8$seeds_tf,list("D"=out.scenario.8$D,
                                                                    "Tp_d"=out.scenario.8$Tp_dd,
                                                                    "Co_r"=out.scenario.8$Co_r,
                                                                    "C_prim"=out.scenario.8$C_prim_f),mean))

out.scenario.mean.8 = data.frame(out.scenario.mean.8,out.scenario.mean.tf.8[,5])
names(out.scenario.mean.8)[5:6] = c("mu","tf")
out.scenario.mean.8$rel = (out.scenario.mean.8$mu-out.scenario.mean.8$tf)/(out.scenario.mean.8$tf)*100
out.scenario.mean.8 = out.scenario.mean.8[out.scenario.mean.8$D==8,]

dat= data.frame(y=0)
ggplot() +
  geom_line(aes(x=Tp_d,y=rel,colour=factor(C_prim)),data=out.scenario.mean.8) +
  facet_grid(Co_r~.,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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



######################### Analysis 6
# vary feeding rate

# set parameters
fPrim = c(0)
C_def = 8
D = c(1,8)
Ths = 11
The = c(180)
Tps = 3
C_prim =  seq(from=0,to=1,length.out=15)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = seq(from=0,to=35,length.out=36)*1e-08
Thd = 6
Tpd = 8
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Ea","Fa")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
#a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s,(a$Tp_d+a$Tp_s)) - a$Tp_s
#a = unique(a[,-15])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Ft = 80
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_varying_Cor.csv",sep=""))

# weights
counts = unique(pars$Co_r/1e-08)
vec = seq(0.01,0.99,length.out=20)
dnbinom(counts,mu=7.67*exp(-(log(exp(3.39))-3.39)^2/(2*0.6167^2)),size=1.224)
dnbinom(38,mu=7.67*exp(-(log(exp(3.39))-3.39)^2/(2*0.6167^2)),size=1.224)
# read

## analyse mean feeding rate..
out.scenario = read.csv("20160801_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_varying_Cor.csv")
out.scenario$rel = (out.scenario$seeds_mu-out.scenario$seeds_tf)/(out.scenario$seeds_tf)*100



ggplot() +
  geom_line(aes(x=Co_r,y=rel,colour=factor(C_prim_f)),data=out.scenario) +
  facet_grid(D+responsepat~alpha,scales="free") +
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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
        #scale_x_log10()+
  ylab("Relative change (%)")

out.scenario$weights = dnbinom(out.scenario$Co_r/1e-08,mu=7.67*exp(-(log(exp(3.39))-3.39)^2/(2*0.6167^2)),size=1.224)

out = melt(lapply(split(out.scenario, list(out.scenario$C_prim_f,out.scenario$D,out.scenario$alpha,out.scenario$responsepat)), function(z) weighted.mean(z$rel, z$weights)) )

# vary probability of attack
out.scenario = read.csv("20160801_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_varying_Cor.csv")
out.scenario0 = out.scenario[out.scenario$Co_r ==0,]
out.scenario = out.scenario[out.scenario$Co_r >0,]
out.scenario$fac = paste(out.scenario$D,out.scenario$C_prim_f,out.scenario$responsepat,out.scenario$alpha)
out.scenario$fac1 = paste(out.scenario$Co_r,out.scenario$D,out.scenario$C_prim_f,out.scenario$responsepat,out.scenario$alpha)

out.scenario0$fac = paste(out.scenario0$D,out.scenario0$C_prim_f,out.scenario0$responsepat,out.scenario0$alpha)

out.scenario$mu0  = out.scenario0[match(out.scenario$fac,out.scenario0$fac),"seeds_mu"]
out.scenario$tf0  = out.scenario0[match(out.scenario$fac,out.scenario0$fac),"seeds_tf"]

out.new = expand.grid(Co_r = unique(out.scenario$Co_r),D=unique(out.scenario$D),C_prim = unique(out.scenario$C_prim_f),responsepat = unique(out.scenario$responsepat),alpha = unique(out.scenario$alpha),prob=seq(from=0,to=1,length.out=20))
out.new$fac = paste(out.new$Co_r,out.new$D,out.new$C_prim,out.new$responsepat,out.new$alpha)

out.new[,(ncol(out.new)+1):(ncol(out.new)+4)] = out.scenario[match(out.new$fac,out.scenario$fac1),c("seeds_tf","seeds_mu","mu0","tf0")]
out.new$rel = NA

for (i in 1:nrow(out.new)){
  out.new$rel[i] = 100*(out.new$prob[i]*(out.new[i,"seeds_mu"]-out.new[i,"seeds_tf"]) + (1-out.new$prob[i])*(out.new[i,"mu0"]-out.new[i,"tf0"]))/(out.new$prob[i]*out.new[i,"seeds_tf"] + (1-out.new$prob[i])*out.new[i,"tf0"])
}

out.new.sub = out.new[out.new$alpha ==1 & out.new$D ==8,]
ggplot()+
  geom_line(aes(x=prob,y=rel,colour=factor(Co_r)),data=out.new.sub)+
  facet_grid(C_prim+D~ alpha+responsepat)+
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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


out.new.sub$weights = dnbinom(out.new.sub$Co_r/1e-08,mu=7.67*exp(-(log(exp(3.39))-3.39)^2/(2*0.6167^2)),size=1.224)
out.new.sub$weights =out.new.sub$weights/sum(dnbinom(c(1:35),mu=7.67*exp(-(log(exp(3.39))-3.39)^2/(2*0.6167^2)),size=1.224))
out.new.sub.mean =   melt(lapply(split(out.new.sub, list(out.new.sub$prob,out.new.sub$C_prim,out.new.sub$D,out.new.sub$alpha,out.new.sub$responsepat)), function(z) weighted.mean(z$rel, z$weights)) )

out.names = t(sapply(split(out.new.sub, list(out.new.sub$prob,out.new.sub$C_prim,out.new.sub$D,out.new.sub$alpha,out.new.sub$responsepat)),function(x){c(x$prob[1],x$D[1],x$C_prim[1],as.character(x$responsepat[1]),x$alpha[1])}))
out.new.sub1 = data.frame(out.new.sub.mean[,1],out.names[,1:ncol(out.names)])
names(out.new.sub1) = c("rel","prob","D","C_prim","responsepat","alpha")
out.new.sub1$prob = as.numeric(out.new.sub1$prob)/20
out.new.sub1 = out.new.sub1[out.new.sub1$alpha ==1 & out.new.sub1$D ==8,]
library(scales)
ggplot()+
  geom_line(aes(x=prob,y=rel,colour=factor(C_prim)),data=out.new.sub1)+
  facet_grid(D~ alpha+responsepat)+
  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
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
  scale_fill_discrete(labels = scales::unit_format("k", 2))+
  ylab("Relative change (%)")


#### simulations 7 Costs of priming when no attack occurs for different durations of priming....

# costs and benefits of priming when different priming durations

# set parameters
fPrim = c(0)
C_def = 8
D = c(8)
Ths = 50# seq(from=1,to=16,1)
The = c(180)
Tps = c(3,10)
C_prim =  seq(0,1,length.out = 50)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = 0#c(6*1e-08,22*1e-08)
Thd = 6
Tpd = seq(from=1,to=16,1)
alpha = c(1)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Th_s = a$Th_s + a$Tp_s
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
#a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s,(a$Tp_d+a$Tp_s)) - a$Tp_s
#a = unique(a[,-15])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Ft = 80
# set herbivory and defense to zero.
pars[,c("Def.f","Def.r","H.f","H.r")] = FALSE
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying_optimal_Tpd_Cor0_Cprim0a.csv",sep=""))



###### Time dynamics and costs of priming

# set parameters
fPrim = c(0)
C_def = 8
D = c(8)
Ths = c(10)
The = c(180)
Tps = c(-7)
C_prim =  c(0.1)
Co_r = c(0,1*1e-8,11*1e-08,22*1e-08)
Thd = 6
Tpd = 100
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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

pars[pars$Co_r ==0,c("H.f","H.r","Prim.f","Prim.r","Def.f","Def.r")] = FALSE

pars$Ft = 80
# do simulations
pars$feeding = "lognorm"
results1 = find_endpoints_defense(pars,Ft=F)
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_Dynamics",sep=""))

#pars = pars[pars$C_prim_f ==0.1,]

out = NA
for (i in 1:nrow(pars)){
  temp = read.csv(paste("Output_defense_",pars$no[i],"_Defense_linear_8_delay20161220.csv",sep=""))
  temp$Co_r = pars$Co_r[i]
  temp$alpha = pars$alpha[i]
  temp$C_prim_f = pars$C_prim_f[i]
  temp$Tp_s = pars$Tp_s[i]
  temp$Th_s = pars$Th_s[i]
  out= rbind(out,temp)
}

out = na.omit(out)
out[out$X ==1,] = 0
out$t = out$X/24
out = out[!out$X==0,]
out$factors = as.factor(paste(out$alpha,out$Co_r,out$C_prim_f,sep=" - "))
out.1=out[out$alpha ==1 & out$Co_r %in% c(0,2.2e-07),]
out.1 = melt(out.1[,c("LAI_f","LAI_mu_f","factors","t")],id.var=c("t","factors"))
names(out.1) = c("t","factors","Strategy","value")
levels(out.1$Strategy) = c("LAI naive","LAI primed")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols= gg_color_hue(8)
cols1 = cols[c(1,4)]
p1 = ggplot()+
  geom_line(aes(x=t,y=value,colour=factors,linetype=Strategy),data=out.1)+
  geom_vline(xintercept=3,linetype=2,size=0.2)+
  geom_vline(xintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=16,linetype=2,size=0.2)+
  geom_text(aes(x=3, label="T[ps]", y=0.25), colour="blue", angle=0,parse=T)+
  geom_text(aes(x=10, label="T[hs]", y=0.25), colour="red", angle=0,parse=T)+
  geom_text(aes(x=16, label="T[hs]+T[hd]", y=0.25),hjust=0.2, colour="red", angle=0,parse=T)+
  ylab("LAI of focal plant (-)")+
  xlab("time (days)")+
  scale_colour_manual("Alpha - Feeding Rate - C_prim",values=cols1)+
  scale_y_continuous()+
  theme(legend.position="right",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  guides(colour= guide_legend(title.position="top"))
#dev.off()

#out = out[out$C_prim_f !=0,]

#out$rel = (out$LAI_mu_f-out$LAI_f)/out$LAI_f
out$rel = (out$Pn_mu_f-out$Pn_f)/out$Pn_f
out$t = out$X/24


out$factors = paste(out$alpha,out$Co_r,out$C_prim_f,sep=" - ")

#out = out[!out$fac =="0 - 0 - 0",]
#pdf("Benefits_costs_priming.pdf")
p2= ggplot()+
  geom_line(aes(x=t,y=rel,colour=factor(factors)),data=out)+
  geom_vline(xintercept=3,linetype=2,size=0.2)+
  geom_vline(xintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=16,linetype=2,size=0.2)+
  geom_text(aes(x=3, label='T[ps]', y=0.1), colour="blue", angle=0,parse=T)+
  geom_text(aes(x=10, label="T[h_s]", y=0.1), colour="red", angle=0,parse=T)+
  geom_text(aes(x=16, label="T[hs]+T[hd]", y=0.1),hjust=0.2, colour="red", angle=0,parse=T)+
  ylab("Relative difference in LAI (-)")+
  xlab("time (days)")+
  scale_y_continuous(limits=c(-0.05,0.1))+
  scale_colour_manual("Alpha - Feeding Rate - C_prim",values=cols)+
  theme(legend.position="right",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  guides(colour= guide_legend(title.position="top"))
#dev.off()

out.p4 = melt(out[,c("LAI_f","LAI_mu_f","factors","t")],id.var=c("t","factors"))
names(out.p4) = c("t","factors","Strategy","value")
levels(out.p4$Strategy) = c("LAI naive","LAI primed")

p4 = ggplot()+
  geom_line(aes(x=t,y=value,colour=factors,linetype=Strategy),data=out.p4)+
  geom_vline(xintercept=3,linetype=2,size=0.2)+
  geom_vline(xintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=16,linetype=2,size=0.2)+
  geom_text(aes(x=3, label="T[ps]", y=0.25), colour="blue", angle=0,parse=T)+
  geom_text(aes(x=10, label="T[hs]", y=0.25), colour="red", angle=0,parse=T)+
  geom_text(aes(x=16, label="T[hs]+T[hd]", y=0.25),hjust=0.2, colour="red", angle=0,parse=T)+
  ylab("LAI of focal plant (-)")+
  xlab("time (days)")+
  scale_colour_manual("Alpha - Feeding Rate - C_prim",values=cols)+
  theme(legend.position="bottom",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  guides(colour= guide_legend(title.position="top"))

# extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p4)


pdf("Combined_dynamics.pdf",width=12,height=8)
p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 2))

dev.off()

### calculate ecological costs of priming during priming event
out = na.omit(out)
out$t = out$X/24
out = out[out$Co_r==2.2e-07,]
out = out[out$C_prim_f ==0.1 & out$alpha==1,]
#
out$MetPrimCosts = out$LAI_mu_f * out$C_prim_f * 0.125 * 8 
out$D.Pn = out$Pn_f - out$Pn_mu_f
# photosynthesis of non-priming strategy minus photosynthesis of priming strategy minus priming costs
out$EcolCosts = out$Pn_f - out$Pn_mu_f -out$MetPrimCosts
priming.temp = out[out$t >= out$Tp_s & out$t < out$Th_s,]  

eco = tapply(priming.temp$EcolCosts,list(priming.temp$alpha,priming.temp$Co_r),sum)
pri = tapply(priming.temp$MetPrimCosts,list(priming.temp$alpha,priming.temp$Co_r),sum)
eco/(pri+eco)

priming.temp$`Direct Costs` = priming.temp$MetPrimCosts/priming.temp$LAI_mu_f
priming.temp$`Indirect Costs` = priming.temp$EcolCosts/priming.temp$LAI_mu_f
priming.temp$alpha = factor(priming.temp$alpha)
priming.temp.2 = melt(priming.temp[,c("t","alpha","Indirect Costs","Direct Costs")],id.var=1:2)
priming.temp.2$t = priming.temp.2$t-3
names(priming.temp.2)[3] ="Costs"



pdf("Ecological costs priming.pdf",width=8,height=6)
ggplot()+
  geom_line(aes(x=t,y=value,colour=Costs,linetype=alpha),data=priming.temp.2)+
  labs(y=expression(paste("Direct/Indirect Costs Priming per unit LAI (",mu,"mol/m"^{2},"s)")))+
  theme(legend.position="right",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  guides(colour= guide_legend(title.position="top"))+
  scale_x_continuous("time since priming (days)")
dev.off()

summary(lm(priming.temp$`Ecological Costs` ~priming.temp$t*priming.temp$alpha))

##################################################
# explore the effect of plant age with herbivory being a constant fraction of LAI

# set parameters
fPrim = c(0)
C_def = 8
D = c(1,2,8)
Ths = c(3,10,20,30)
The = c(180)
Tps = c(-2,-8,-12)
C_prim =  c(0,0.1,0.25,0.5,0.75,1)
Co_r = c(0.6e-07,1.1e-07,2.2e-07)
Thd = 6
Tpd = 100
alpha = c(1,1.25)
feeding = c("prop")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
# temp = pars[pars$fPrim_f ==0 & pars$Co_r ==1.1e-07,]
# pars$ref_Co_r = pars$Co_r
# temp$ref_Co_r = temp$Co_r
# temp$Co_r = 0
# temp[,c("Def.f","Def.r","H.f","H.r")] = FALSE
# temp$no = seq(from=(max(pars$no)+1),to=(nrow(temp)+max(pars$no)))
# pars = rbind(pars,temp)
# optimal flowering time without herbivory
Opt.flower = data.frame(D=c(1,2,6,8),Ft=c(90.58333333,84.875,75.83333333,73.45833333))
pars$Ft = Opt.flower[match(pars$D,Opt.flower$D),2]
pars$Ft = 80
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft.csv",sep=""))

out.scenario = read.csv("20160726_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft.csv")
out.scenario$Tp_ss = out.scenario$Tp_s - out.scenario$Th_s

# relative benefit of priming compared to non-priming
out.scenario$rel = -out.scenario$seeds_tfm/(out.scenario$seeds_tf)*100


##### sensitivity analysis for different parameters.
### The following parameters are changed while keeping the others at a specified level:
#C_prim (done see analysis)
#alpha
#feedingrate (done)
#herbivore onset 
#delay in defense of non-primed plants 
#time between priming event and trigger event (done see read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying.csv"))


######## C_prim

    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = c(10)#c(3,10,20,30)
    The = c(180)
    Tps = c(-8)#c(-2,-8,-12)
    C_prim =  seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    #Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
    Co_r = c(8*1e-08)#c(22*1e-08)
    Thd = 6
    Tpd = 100
    alpha = c(1)#c(1,1.1,1.25)
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])
    
    
    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    
    # if Co_r is 0, Def.f and Def.r FALSE
    pars$Ft = 80
    # do simulations
    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    #write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cprim.csv",sep=""))

    result = read.csv("20161220_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cprim.csv")
    result$rel = -result$seeds_tfm/(result$seeds_tf)*100
    result$Density = factor(result$D)
p1= ggplot() +
      geom_line(aes(x=C_prim_f,y=rel,colour=Density),data=result) +
    #  geom_hline(aes(yintercept=y),data=dat) +
      theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
            axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
            panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
            size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
            size = 0.5), strip.background = element_rect(fill = "grey80", 
            colour = "grey50"), legend.text=element_text(size=12),
            strip.text.x = element_text(size=12),
            strip.text.y = element_text(size=12)) +
      ylab("")+
      #ylab("Relative difference in seed production (%)")+
      xlab("Priming level (fPrim [-])")

 
 ######## alpha
    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = c(10)#c(3,10,20,30)
    The = c(180)
    Tps = c(-8)#c(-2,-8,-12)
    C_prim =  0.1#seq(0,0.5,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    #Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
    Co_r = c(8*1e-08)#c(22*1e-08)
    Thd = 6
    Tpd = 100
    alpha = seq(1,1.5,length.out=20)#c(1,1.1,1.25)
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])
    
    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    
    # if Co_r is 0, Def.f and Def.r FALSE
    pars$Ft = 80

    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    #write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
    
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_alpha.csv",sep=""))
    
    
     result = read.csv("20170501_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_alpha.csv")
     result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$D = factor(result$D)
 p2= ggplot() +
       geom_line(aes(x=alpha,y=rel,colour=D),data=result) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
             axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
             panel.background = element_rect(fill = "white",
             colour = NA), panel.border = element_rect(fill = NA,
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90",
             size = 0.2), panel.grid.minor = element_line(colour = "grey98",
             size = 0.5), strip.background = element_rect(fill = "grey80",
             colour = "grey50"), legend.text=element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12)) +
     #  ylab("Relative difference in seed production (%)")+
       ylab("")+
       scale_y_continuous(limits=c(0,100))+ 
       xlab("Competition coefficient (alpha [-])")
# 
    
######## feedingrate
    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = c(10)#c(3,10,20,30)
    The = c(180)
    Tps = c(-8)#c(-2,-8,-12)
    C_prim =  0.1#seq(0,1,length.out=20)#0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    Co_r = seq(from=0,to=35,length.out=36)*1e-08
    Thd = 6
    Tpd = 100
    alpha = 1
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])
    
    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    pars$Ft = 80
    
    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cor.csv",sep=""))

     result = read.csv("20161220_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cor.csv")
   result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$D = factor(result$D)
 p3=  ggplot() +
       geom_line(aes(x=Co_r,y=rel,colour=D),data=result) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
             axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12)) +
       ylab("Relative difference in seed production (%)")+
       xlab(expression(paste("Feeding rate (Cor [m"^{2},"s])")))
    

###### herbivore onset
    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = seq(9,30,length.out=22)
    The = c(180)
    Tps = c(-8)#c(-2,-8,-12)
    C_prim =  0.1#seq(0,1,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    Co_r = 8e-08#c(22*1e-08)
    Thd = 6
    Tpd = 100
    alpha = 1
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])

    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    pars$Ft = 80
    
    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    
    
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Ths.csv",sep=""))

     result = read.csv("20170501_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Ths.csv")
     result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$Density = factor(result$D)
 p4= ggplot() +
       geom_line(aes(x=Th_s,y=rel,colour=Density),data=result) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
             axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5),  strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12)) +
       ylab("")+
#   #    ylab("Relative difference in seed production (%)")+
       xlab(expression(paste("Start of trigger event (",T[hs],",[days])")))
    

######## delay in defense of non-primed plants

    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = 8
    The = c(180)
    Tps = c(-8)#c(-2,-8,-12)
    C_prim =  0.1#seq(0,0.5,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    Co_r = 8*1e-08#c(22*1e-08)
    Thd = seq(0,14,length.out = 15*4)
    Tpd = 100
    alpha = 1
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])
    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    pars$Ft = 80
    
    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Thd.csv",sep=""))

     result = read.csv("20170501_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Thd.csv")
     result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$D = factor(result$D)
 p5= ggplot() +
       geom_line(aes(x=Th_d,y=rel,colour=D),data=result) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
             axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5),  strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12)) +
             xlim(c(0,10))+
      # ylab("Relative difference in seed production (%)")+
       ylab("")+
       xlab(expression(paste("Delay in defense non primed plants (",T[hd]," [days]")))

######## time between priming and trigger event

    fPrim = c(0)
    C_def = 8
    D = c(1,2,6,8)
    Ths = 8
    The = c(180)
    Tps = seq(-7,-1,by=0.5)
    C_prim =  0.1#seq(0,1,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
    Co_r = 8*1e-8#c(22*1e-08)
    Thd = 6
    Tpd = 100
    alpha = 1
    feeding = c("lognorm")
    responsepat = c("Ea")
    # make parameter matrix
    a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
    a$Tp_ss = a$Tp_s
    a$Tp_dd = a$Tp_d
    a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
    a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
    a = unique(a[,-14])
    
    # add to full parameter dataframe
    comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
    pars$Th_d = a$Th_d
    pars$Tco_d = 12
    pars$shape = "delay"
    pars$alpha = a$alpha
    pars$feeding = a$feeding
    pars$responsepat = a$responsepat
    pars$Ft = 80
    
    results1 = find_endpoints_defense(pars,Ft=F)
    results1 = data.frame(pars,results1)
    results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
    # store results
    
    write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Tps.csv",sep=""))
    
     result = read.csv("20170501_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Tps.csv")
     result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$D = factor(result$D)
 p6=    ggplot() +
       geom_line(aes(x=-Tp_s,y=rel,colour=D),data=result) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
             axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12)) +
       ylab("")+
     #  ylab("Relative difference in seed production (%)")+
       xlab("Time priming and trigger event (days)")

# extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p4)


plots1 <- list(p1 + theme(legend.position="none"),
               p2 + theme(legend.position="none"))
plots2 <- list(p3 + theme(legend.position="none"),
               p4 + theme(legend.position=c(0.8,0.5),
                          legend.key=element_rect(fill = "white", colour = NA,size=12),
                          legend.title= element_text(size=12),
                          legend.background =element_rect(fill = "white", colour = NA)))
plots3 <- list(p5 + theme(legend.position="none"),
               p6 + theme(legend.position="none"))
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
grobs3 = lapply(plots3, ggplotGrob)

g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(cbind, c(grobs3, size="first"))
pdf("Ecological factors_poster.pdf")
g4 = do.call(rbind, c(list(g1,g2,g3), size="first")) #combine g1 and g2 into a list
grid.draw(g4)
dev.off()

plots <- list(p1,p2,p3,p4)
grobs = lapply(plots, ggplotGrob)
g = do.call(rbind, c(grobs, size="first"))

g$widths = do.call(unit.pmax, lapply(grobs, "[[", "widths"))
grid.newpage()
grid.draw(g)


grid.arrange(grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2),
                ggplotGrob(p3), ggplotGrob(p4),
                ggplotGrob(p5), ggplotGrob(p6),
                size="first")),ncol=2)


# analyse the effect of ecological factors on the optimal level of priming 

#C_prim (done see analysis)

# Not applicable

#alpha

     result = read.csv("20161222_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_alpha.csv")
     result$rel = -result$seeds_tfm/(result$seeds_tf)*100
     result$Density = factor(result$D)
     index= melt(tapply(result$rel,list(result$alpha,result$D),which.max))
     names(index) = c("alpha","Density","value")
     index$Density =factor(index$Density)
     C_prim  = unique(result$C_prim_f)
     index$fPrim = C_prim[index$value]
     
 p2= ggplot() +
       geom_line(aes(x=alpha,y=fPrim,colour=Density),data=index) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=9), axis.text.y=element_text(size=9),
             axis.title.y=element_text(size=9),axis.title.x=element_text(size=9),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=9),
             strip.text.x = element_text(size=9),
             strip.text.y = element_text(size=9)) +
       ylab("")+
   #    ylab("fPrim")+
       xlab("alpha")


#feedingrate (done)
 
 result = read.csv("20161221_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cor.csv")
 result$rel = -result$seeds_tfm/(result$seeds_tf)*100
 result$Density = factor(result$D)
 index= melt(tapply(result$rel,list(result$Co_r,result$D),which.max))
 names(index) = c("Co_r","Density","value")
 index$Density =factor(index$Density)
 C_prim  = unique(result$C_prim_f)
 index$fPrim = C_prim[index$value]
 
p3= ggplot() +
   geom_line(aes(x=Co_r,y=fPrim,colour=Density),data=index) +
   #  geom_hline(aes(yintercept=y),data=dat) +
   theme(axis.text.x=element_text(angle=0,size=9), axis.text.y=element_text(size=9),
         axis.title.y=element_text(size=9),axis.title.x=element_text(size=9),
         panel.background = element_rect(fill = "white", 
         colour = NA), panel.border = element_rect(fill = NA, 
         colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
         size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
         size = 0.5), strip.background = element_rect(fill = "grey80", 
         colour = "grey50"), strip.background = element_rect(fill = "grey80", 
         colour = "grey50"), legend.text=element_text(size=9),
         strip.text.x = element_text(size=9),
         strip.text.y = element_text(size=9)) +
   ylab("")+
   #    ylab("fPrim")+
   xlab(expression(paste("Feeding rate (Cor [",mu,"mol/m"^{2},"s])")))
 
#herbivore onset 

    result = read.csv("20161222_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Ths.csv")
    result$rel = -result$seeds_tfm/(result$seeds_tf)*100
    result$Density = factor(result$D)
    index= melt(tapply(result$rel,list(result$Th_s,result$D),which.max))
    names(index) = c("Ths","Density","value")
    index$Density =factor(index$Density)
    C_prim  = unique(result$C_prim_f)
    index$fPrim = C_prim[index$value]

p4= ggplot() +
      geom_line(aes(x=Ths,y=fPrim,colour=Density),data=index) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=9), axis.text.y=element_text(size=9),
             axis.title.y=element_text(size=9),axis.title.x=element_text(size=9),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=9),
             strip.text.x = element_text(size=9),
             strip.text.y = element_text(size=9)) +
       ylab("")+
   #    ylab("Relative difference in seed production (%)")+
       xlab(expression(paste("Start of trigger event (",T[hs],",[days])")))


#delay in defense of non-primed plants 

    result = read.csv("20161222_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Thd.csv")
    result$rel = -result$seeds_tfm/(result$seeds_tf)*100
    result$Density = factor(result$D)
    index= melt(tapply(result$rel,list(result$Th_d,result$D),which.max))
    names(index) = c("Thd","Density","value")
    index$Density =factor(index$Density)
    C_prim  = unique(result$C_prim_f)
    index$fPrim = C_prim[index$value]

p5= ggplot() +
  geom_line(aes(x=Thd,y=fPrim,colour=Density),data=index) +
  #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=9), axis.text.y=element_text(size=9),
             axis.title.y=element_text(size=9),axis.title.x=element_text(size=9),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=9),
             strip.text.x = element_text(size=9),
             strip.text.y = element_text(size=9)) +
             xlim(c(0,10))+
   #    ylab("Relative difference in seed production (%)")+
       ylab("")+
       xlab(expression(paste("Delay in defense non primed plants (",T[hd]," [days]")))

#time between priming event and trigger event (done see read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying.csv"))

    result = read.csv("20161222_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Tps.csv")
    result$rel = -result$seeds_tfm/(result$seeds_tf)*100
    result$Density = factor(result$D)
    index= melt(tapply(result$rel,list(result$Tp_s,result$D),which.max))
    names(index) = c("Tp_ss","Density","value")
    index$Density =factor(index$Density)
    C_prim  = unique(result$C_prim_f)
    index$fPrim = C_prim[index$value]

p6=    ggplot() +
       geom_line(aes(x=-Tp_ss,y=fPrim,colour=Density),data=index) +
       #  geom_hline(aes(yintercept=y),data=dat) +
       theme(axis.text.x=element_text(angle=0,size=9), axis.text.y=element_text(size=9),
             axis.title.y=element_text(size=9),axis.title.x=element_text(size=9),
             panel.background = element_rect(fill = "white", 
             colour = NA), panel.border = element_rect(fill = NA, 
             colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
             size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
             size = 0.5), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), strip.background = element_rect(fill = "grey80", 
             colour = "grey50"), legend.text=element_text(size=9),
             strip.text.x = element_text(size=9),
             strip.text.y = element_text(size=9)) +
       ylab("")+
     #  ylab("Relative difference in seed production (%)")+
       xlab("Time priming and trigger event (days)")

# function to extract a column from a dataframe with name tat is stored in mod
fun.name = function(x,mod){
  a = x[,eval(mod)]
  return(a)
}



par.test$alpha
fun.name(par.test,"alpha")
a = "alpha"
te = paste("par.test$",a,sep="")
eval(te)
get(te)
par.test[,eval(a)] = 1
par.test[,"alpha"]
assign(eval(te),00.11)
a <- 1:4
assign("a[1]", 2,envir=.GlobalEnv)
a[1] == 2          # FALSE
get("a[1]") == 2   # TRUE
eval(parse(te))

par.test= pars[1,]
parvalue=0.1

par.test$alpha=1.5
par.test$D=8
optimal = function(par,other.pars,rel=T){
    other.pars$C_prim_f = par
    result = find_endpoints_defense(other.pars,Ft=F)  
    if (rel){
      out = result$seeds_tfm/(result$seeds_tf)*100  
    } else {
      out = result$seeds_tfm
    }
    
    return(out)
}

optimal(par=10.5,par.test,rel=F)

optimal.optimize = function(par,other.pars){
  other.pars$C_prim_f = par
  result = find_endpoints_defense(other.pars,Ft=F)  
  # this gives the negative fitness for minimization
#  result$rel = result$seeds_tfm/(result$seeds_tf)*100
#  return(result$rel)
  return(result)
}


par = 0.1
out = optim(par,optimal,other.pars=par.test,method="Brent",lower=0,upper=1,control=list(trace=1,REPORT=1))

a = "alpha"
par.test[,eval(a)] = 1.5

##############################

#### intialize paras

fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(10)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  0.1# seq(0,0.5,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
#Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Co_r = c(22*1e-08)
Thd = 6
Tpd = 100
alpha = 1
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "delay"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80
pars$C_def_max = pars$C_def

par.test = pars[1,]
par.test$Tp_ss = -4

alpha = data.frame(par="alpha",value=seq(1,1.5,length.out=20))
Co_r = data.frame(par="Co_r",value=seq(from=0,to=35,length.out=36)*1e-08)
Ths = data.frame(par="Th_s",value=seq(5,25,by=1))
Thd = data.frame(par="Th_d",value=seq(0,14,length.out = 29))
Tps = data.frame(par="Tp_s",value=seq(1,15,by=0.5))
Ft = data.frame(par="Th_e",value=seq(60,100,length.out=21))
C_def = data.frame(par="C_def",value=seq(0.5,8,length.out=20))
#par.var = rbind.data.frame(Tps,Thd,Ths,Co_r,alpha)
par.var = rbind.data.frame(C_def)
#D=c(1,2,6,8)
#name = expand.grid(D,par.var[,1])
#par.var.e = expand.grid(D=D,value=par.var[,2])
#responsepat = c("Ea","Fa")
responsepat = c(Co_r = c(8e-08))
name = expand.grid(responsepat,par.var[,1],Th_d =c(6,10))
par.var.e = expand.grid(D=responsepat,value=par.var[,2],Th_d =c(6,10))
par.var.e$parname = as.character(name[,2])
par.var.e$rel = NA
par.var.e$fPrim = NA
par.var.e$converge = NA
par.var.e$ir = NA
par.var.e$seeds_tf = NA
#par.var.e = par.var.e[par.var.e$D==8 & par.var.e$parname == "Tp_s",]
# set parameters
par.test$alpha = 1.0
par.test$Tp_d = 20
par.test$Th_s = 15
par.test$D = 8
par.test$shape = "linear"

write.table(1,"par.var.e.responspat_Dvar.csv",sep=",")
for (i in 1:nrow(par.var.e)){#nrow(par.var.e)){
  if (ncol(par.var.e)<8) print("adjust size of par.var.e")
  par.temp = par.test
  par.temp[,eval(par.var.e[i,4])] = par.var.e[i,2]
  par.temp$C_def_f = par.temp$C_def_r = par.temp$C_def_max = par.temp$C_def
  # if parameter to be changed is Th_s then also the day when priming starts changes. 
  if (eval(par.var.e[i,4]) == "Th_s"){
    par.temp$Tp_s = par.temp$Th_s + par.test$Tp_ss
  }
  if (eval(par.var.e[i,4]) == "Tp_s"){
    par.temp$Th_s = 15
  }
  # set duration of priming such that it is maximal the time between tps and ths
  par.temp$Tp_d = min(par.temp$Tp_s+par.temp$Tp_d,par.temp$Th_s) - par.temp$Tp_s
  # set plant density
  #par.temp[,"D"] = par.var.e[i,1]
  # set responspattern
  #par.temp[,"responsepat"] = par.var.e[i,1]
  # set feeding rate
  par.temp[,"Co_r"] = par.var.e[i,1]
  # set th_d 
  par.temp[,"Th_d"] = par.var.e[i,3]
  par = c(0.5)
  #out1 = optim(par,optimal,lower=c(0),upper=c(1),other.pars=par.temp,method="Brent")
  out = optimize(optimal,lower=c(0),upper=c(1),other.pars=par.temp,maximum=F,rel=T)
  par.var.e$rel[i] = out$objective
  par.var.e$fPrim[i] = out$minimum
  par.var.e$converge[i] = "tt" #out$convergence
  # compare with induced resistance (C_prim =1)
  par.temp$C_prim_f = par.temp$C_prim_r = 1
  out.ir = find_endpoints_defense(par.temp,Ft=F)  
  par.var.e$ir[i] =out.ir$seeds_tfm/(out.ir$seeds_tf)*100  
  par.var.e$seeds_tf[i] = out.ir$seeds_tf
  print(data.frame(Density= par.var.e[i,1],par.temp))
  write.table(par.var.e[i,],"par.var.e.responspat_Dvar.csv",sep=",",append=T,col.names=F)
  print(paste("+++++++++++++++++++++++++",i,"++++++++++++++++++++++++++++++++++++++++",sep="  "))
} 
write.csv(par.var.e,paste(format(Sys.time(), "%Y%m%d"),"_par.var.e.par.var.e.responsepat_Dvar_3.csv",sep=""))


par.var.e = read.csv("20170126_par.var.e.par.var.e.responsepat_D8.csv")
ggplot(data=par.var.e) +
  geom_line(aes(x=value,y=fPrim,colour=factor(D))) +
              facet_grid(.~ parname,scales="free")

par.var.e$Response = factor(par.var.e$D)
par.var.e.alpha = par.var.e[par.var.e$parname=="alpha",]
par.var.e.cor = par.var.e[par.var.e$parname=="Co_r",]
par.var.e.thd = par.var.e[par.var.e$parname=="Th_d",]
par.var.e.ths = par.var.e[par.var.e$parname=="Th_s",]
par.var.e.tps = par.var.e[par.var.e$parname=="Tp_s",]
par.var.e.tps$value = 15-par.var.e.tps$value

p1= ggplot() + 
  geom_line(aes(x=value,y=fPrim,colour=Response),data=par.var.e.alpha)+
         theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
               axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
               panel.background = element_rect(fill = "white", 
               colour = NA), panel.border = element_rect(fill = NA, 
               colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
               size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
               size = 0.5), strip.background = element_rect(fill = "grey80", 
               colour = "grey50"), strip.background = element_rect(fill = "grey80", 
               colour = "grey50"), legend.text=element_text(size=12),
               strip.text.x = element_text(size=12),
               strip.text.y = element_text(size=12)) +
       ylab("")+
#       ylab("Priming level (fPrim [-])")+
       xlab("Competition coefficient (alpha [-])")

  
p2= ggplot() + 
  geom_line(aes(x=value,y=fPrim,colour=Response),data=par.var.e.cor)+
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
       ylab("")+
#  ylab("Priming level (fPrim [-])")+
  xlab(expression(paste("Feeding rate (",C[or],"[m"^{2},"s])")))

p3= ggplot() + 
  geom_line(aes(x=value,y=fPrim,colour=Response),data=par.var.e.thd)+
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
  #       ylab("")+
  ylab("Priming level (fPrim [-])")+
  xlab(expression(paste("Delay in defense non primed plants (",T[hd]," [days]")))

p4= ggplot(aes(x=value,y=fPrim,colour=Response),data=par.var.e.ths) + 
  geom_line(aes(x=value,y=fPrim))+
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
     ylab("")+
#  ylab("Priming level (fPrim [-])")+
  xlab(expression(paste("Start of trigger event (",T[hs],",[days])")))

p5= ggplot() + 
  geom_line(aes(x=value,y=fPrim,colour=Response),data=par.var.e.tps)+
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
       ylab("")+
#  ylab("Priming level (fPrim [-])")+
  xlab(expression(paste("Time between priming and trigger event (",T[ps]," [days]")))

dat=data.frame(x=c(0,10),y=c(0,10))
#eliminates background, gridlines, and chart border
p6= ggplot(data=dat)+
  #plots the points
  geom_point(aes(x=x,y=y),colour="NA") +
    #theme with white background
  theme_bw() +
  
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank(),
    axis.text = element_text(colour="white")
  ) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'white'), 
        axis.ticks.length=unit(0,"cm"))+
  xlab("") + ylab("")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p3)


plots1 <- list(p1 + theme(legend.position=c(0.25,0.5),
                          legend.key=element_rect(fill = "white", colour = NA,size=12),
                          legend.title= element_text(size=12),
                          legend.background =element_rect(fill = "white", colour = NA)),
               p2 + theme(legend.position="none"))
plots2 <- list(p3 + theme(legend.position="none"),
               p4 + theme(legend.position="none"))
plots3 <- list(p5 + theme(legend.position="none"),
               p6+ theme(legend.position="none"))
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
grobs3 = lapply(plots3, ggplotGrob)


g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(cbind, c(grobs3, size="first"))
pdf("Optimal fPrim_Eafa_poster.pdf")
g4 = do.call(rbind, c(list(g1,g2,g3), size="first")) #combine g1 and g2 into a list
grid.draw(g4)
dev.off()
  
plots <- list(p1,p2,p3,p4)
grobs = lapply(plots, ggplotGrob)
g = do.call(rbind, c(grobs, size="first"))

g$widths = do.call(unit.pmax, lapply(grobs, "[[", "widths"))
grid.newpage()
grid.draw(g)

##############################################################
################ revision manuscript ########################
##############################################################

# Ecological factors: Fig. 5

##### sensitivity analysis for different parameters.
### The following parameters are changed while keeping the others at a specified level:
#C_prim (done see analysis)
#alpha
#feedingrate (done)
#herbivore onset 
#delay in defense of non-primed plants 
#time between priming event and trigger event (done see read.csv("20160729_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_fixed_Cor_Tpd_varying.csv"))


######## C_prim

fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(10)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
#Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Co_r = c(8*1e-08)#c(22*1e-08)
Thd = 6
Tpd = 100
alpha = c(1)#c(1,1.1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$C_def_max = pars$C_def

# if Co_r is 0, Def.f and Def.r FALSE
pars$Ft = 80
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cprim.csv",sep=""))

result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cprim.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$Density = factor(result$D)
p1= ggplot() +
  geom_line(aes(x=C_prim_f,y=rel,colour=Density),data=result) +
  #  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
  ylab("")+
  #ylab("Relative difference in seed production (%)")+
  xlab(expression(paste("Priming level (",italic(f)[prim]," [-])",sep="")))


######## alpha
fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(10)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  0.1#seq(0,0.5,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
#Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Co_r = c(8*1e-08)#c(22*1e-08)
Thd = 6
Tpd = 100
alpha = seq(1,1.5,length.out=20)#c(1,1.1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat

# if Co_r is 0, Def.f and Def.r FALSE
pars$Ft = 80

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))

write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_alpha.csv",sep=""))


result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_alpha.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$D = factor(result$D)

p2= ggplot() +
   geom_line(aes(x=alpha,y=rel,colour=D),data=result) +
   #  geom_hline(aes(yintercept=y),data=dat) +
   theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
         panel.background = element_rect(fill = "white",
         colour = NA), panel.border = element_rect(fill = NA,
         colour = "grey50"), panel.grid.major = element_line(colour = "grey90",
         size = 0.2), panel.grid.minor = element_line(colour = "grey98",
         size = 0.5), strip.background = element_rect(fill = "grey80",
         colour = "grey50"), legend.text=element_text(size=12),
         strip.text.x = element_text(size=12),
         strip.text.y = element_text(size=12)) +
   #  ylab("Relative difference in seed production (%)")+
   ylab("")+
   scale_y_continuous(limits=c(-3,0.2))+ 
   xlab(expression(paste("Competition coefficient (",alpha," [-])")))
 

######## feedingrate
fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = c(10)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  0.1#seq(0,1,length.out=20)#0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = seq(from=0,to=35,length.out=36)*1e-08
Thd = 6
Tpd = 100
alpha = 1
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cor.csv",sep=""))

result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Cor.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$D = factor(result$D)
p3=  ggplot() +
 geom_line(aes(x=Co_r,y=rel,colour=D),data=result) +
 #  geom_hline(aes(yintercept=y),data=dat) +
 theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
       axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
       panel.background = element_rect(fill = "white", 
         colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
         size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
         size = 0.5), strip.background = element_rect(fill = "grey80", 
         colour = "grey50"), legend.text=element_text(size=12),
         strip.text.x = element_text(size=12),
         strip.text.y = element_text(size=12)) +
   ylab("Relative difference in seed production (%)")+
   xlab(expression(paste("Feeding rate (",italic(F)," [m"^{2},"s])")))


###### herbivore onset
fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = seq(9,30,length.out=22)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  0.1#seq(0,1,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = 8e-08#c(22*1e-08)
Thd = 6
Tpd = 100
alpha = 1
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results


write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Ths.csv",sep=""))

result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Ths.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$Density = factor(result$D)
p4= ggplot() +
 geom_line(aes(x=Th_s,y=rel,colour=Density),data=result) +
 #  geom_hline(aes(yintercept=y),data=dat) +
 theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
       panel.background = element_rect(fill = "white", 
       colour = NA), panel.border = element_rect(fill = NA, 
       colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
       size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
       size = 0.5),  strip.background = element_rect(fill = "grey80", 
       colour = "grey50"), legend.text=element_text(size=12),
       strip.text.x = element_text(size=12),
       strip.text.y = element_text(size=12)) +
 ylab("")+
 #   #    ylab("Relative difference in seed production (%)")+
 xlab(expression(paste("Start of trigger event (",italic(t)[hs],",[days])")))

######## delay in defense of non-primed plants

fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = 8
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  0.1#seq(0,0.5,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = 8*1e-08#c(22*1e-08)
Thd = seq(0,14,length.out = 15*4)
Tpd = 100
alpha = 1
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])
# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results

write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Thd.csv",sep=""))

result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Thd.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$D = factor(result$D)
p5= ggplot() +
  geom_line(aes(x=Th_d,y=rel,colour=D),data=result) +
  #  geom_hline(aes(yintercept=y),data=dat) +
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white",
        colour = NA), panel.border = element_rect(fill = NA,
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90",
        size = 0.2), panel.grid.minor = element_line(colour = "grey98",
        size = 0.5),  strip.background = element_rect(fill = "grey80",
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
  xlim(c(0,10))+
  ylim(c(-2,25))+
  # ylab("Relative difference in seed production (%)")+
  ylab("")+
  xlab(expression(paste("Delay in defense non primed plants (",italic(t)[hd]," [days]")))

######## time between priming and trigger event

fPrim = c(0)
C_def = 8
D = c(1,2,6,8)
Ths = 8
The = c(180)
Tps = seq(-7,-1,by=0.5)
C_prim =  0.1#seq(0,1,length.out=20) #0.1# seq(0,1,length.out=20)#c(0,0.1,0.25,0.5,0.75,1)
Co_r = 8*1e-8#c(22*1e-08)
Thd = 6
Tpd = 100
alpha = 1
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "sigmoid"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80
pars$C_def_max = pars$C_def

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results

write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Tps.csv",sep=""))

result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_Tps_sigmoid.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$D = factor(result$D)
p6=    ggplot() +
   geom_line(aes(x=-Tp_ss_1,y=rel,colour=D),data=result) +
   #  geom_hline(aes(yintercept=y),data=dat) +
   theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
         panel.background = element_rect(fill = "white", 
         colour = NA), panel.border = element_rect(fill = NA, 
         colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
         size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
         size = 0.5), strip.background = element_rect(fill = "grey80", 
         colour = "grey50"), legend.text=element_text(size=12),
         strip.text.x = element_text(size=12),
         strip.text.y = element_text(size=12)) +
   ylab("")+
   #  ylab("Relative difference in seed production (%)")+
   xlab("Time priming and trigger event (days)")

# extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p4)


plots1 <- list(p1 + theme(legend.position="none"),
               p2 + theme(legend.position="none"))
plots2 <- list(p3 + theme(legend.position="none"),
               p4 + theme(legend.position=c(0.8,0.5),
                          legend.key=element_rect(fill = "white", colour = NA,size=12),
                          legend.title= element_text(size=12),
                          legend.background =element_rect(fill = "white", colour = NA)))
plots3 <- list(p5 + theme(legend.position="none"),
               p6 + theme(legend.position="none"))
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
grobs3 = lapply(plots3, ggplotGrob)

g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(cbind, c(grobs3, size="first"))
pdf("Sigmoidal_Defense_8_1.pdf")
g4 = do.call(rbind, c(list(g1,g2,g3), size="first")) #combine g1 and g2 into a list
grid.draw(g4)
dev.off()


####################################################################################
############## Change costs of defense #############################################

######## C_prim

fPrim = c(0)
C_def = seq(0.5,8,length.out=20)
D = c(8)
Ths = c(10)#c(3,10,20,30)
The = c(180)
Tps = c(-8)#c(-2,-8,-12)
C_prim =  c(0.1,1)# c(0,0.1,0.25,0.5,0.75,1)
#Co_r = c(2*1e-08,6*1e-08,11*1e-08,22*1e-08)
Co_r = c(1e-08,8*1e-08)#c(22*1e-08)
Thd = c(4,6,8,10)
Tpd = 100
alpha = c(1)#c(1,1.1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")
# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])

#a$C_prim = 1/a$C_def


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
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
pars$Th_d = a$Th_d
pars$Tco_d = 12
pars$shape = "linear"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$C_def_max = pars$C_def

# if Co_r is 0, Def.f and Def.r FALSE
pars$Ft = 80
# do simulations
results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_varying thd.csv",sep=""))


result = read.csv("20170720_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_lognormal_fixC_def_Ft_Fig2_panelplot_varying thd.csv")
result$rel = -result$seeds_tfm/(result$seeds_tf)*100
result$D = factor(result$D)
result$Co_r =as.factor(result$Co_r)
result$Co_r = relevel(result$Co_r,ref="8e-08")

result$Th_d = factor(result$Th_d)
result$C_prim_f = as.factor(result$C_prim_f)

pdf("Varying cdef and thd.pdf",width=8,height=6)
#p6=    
  ggplot() +
  geom_line(aes(x=C_def,y=rel,colour=as.factor(Th_d),linetype=C_prim_f),data=result) +
  #  geom_hline(aes(yintercept=y),data=dat) +
    facet_grid(.~ factor(Co_r))+
  theme(axis.text.x=element_text(angle=0,size=12), axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),axis.title.x=element_text(size=12),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
    labs(colour=expression(paste(italic(t)[hd]," [days]")),
         linetype=expression(paste(italic(f)[prim])))+
  #ylab("")+
    ylab("Relative difference in seed production (%)")+
  xlab(expression(paste("Costs of defense (",C[def]," [",mu,"mol/m"^{2},"s]")))
dev.off()

head(result)
temp = result[,c("Co_r","C_prim_f","Th_d","C_def","rel")]
temp1 = temp[temp$C_prim_f==1,]
temp2 = temp[temp$C_prim_f==0.1,]
temp3 = data.frame(temp1[,c(1:4)], y=temp1$rel - temp2$rel)

ggplot(data=temp3)+
  geom_line(aes(x=C_def,y=y,colour=as.factor(Th_d), linetype=as.factor(Co_r)))

####################################################### 
#### Dynamics LAI
# set parameters
fPrim = c(0)
C_def = c(2,4,8)
D = c(8)
Ths = c(10)
The = c(180)
Tps = c(-7)
C_prim =  0.1#c(0.1,0.2)
Co_r = c(1*1e-8,8*1e-08)
Thd = 6
Tpd = 100
alpha = c(1,1.25)
feeding = c("lognorm")
responsepat = c("Ea")

# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
pars = as.data.frame(matrix(NA,nrow=nrow(a),ncol=length(comb[1,])))

names(pars) = names(comb)
pars[] = comb[1,]

#a$C_prim = 1/a$C_def

pars$C_def = a$C_def
pars$C_def_max = pars$C_def
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
pars$Th_d = 8
pars$Tco_d = 6
pars$shape = "linear"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat



pars[pars$Co_r ==0,c("H.f","H.r","Prim.f","Prim.r","Def.f","Def.r")] = FALSE

pars$Ft = 80
# do simulations
pars$feeding = "lognorm"

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_propLAI_fixC_def_Ft.csv",sep=""))


out = NA
for (i in 1:nrow(pars)){
  temp = read.csv(paste("Output_defense_",pars$no[i],"_Defense_linear_8_delay20170712.csv",sep=""))
  temp$Co_r = pars$Co_r[i]
  temp$C_def = pars$C_def[i]
  temp$alpha = pars$alpha[i]
  temp$C_prim_f = pars$C_prim_f[i]
  out= rbind(out,temp)
}
out = na.omit(out)
out[out$X ==1,] = 0
out$t = out$X/24
out = out[!out$X==0,]
out$rel = (out$Pn_mu_f-out$Pn_f)/out$Pn_f

out = out[out$alpha==1,]
out$Co_r = factor(out$Co_r)
out$Co_r=relevel(out$Co_r,ref=2)
out$alpha = factor(out$alpha)
out$C_prim_f = factor(out$C_prim_f)
out = out[out$Co_r %in%  c(1*1e-8,8*1e-08),]
out.group1 = out[out$C_prim_f == 0.1,]
out.group2 = out[out$C_prim_f == 0.25,]
out.group3 = out[out$C_prim_f == 0.5,]
thin =seq(0,80,length.out=11)
out.group1.thin = out.group1[out.group1$t %in% thin,]
out.group2.thin = out.group2[out.group2$t %in% thin,]
out.group3.thin = out.group3[out.group3$t %in% thin,]

points = c("0.1" = 1,"0.2"=2,"0.05"=3)
# alternative figure
#pdf("Priming_costs_constant_effect constant_defense costs variable.pdf")
ggplot()+
  geom_line(aes(x=t,y=rel,colour=Co_r,linetype=factor(C_def)),data=out.group1)+
#  geom_line(aes(x=t,y=rel,colour=Co_r,linetype=factor(C_def)),data=out.group2)+
#  geom_line(aes(x=t,y=rel,colour=Co_r,linetype=factor(C_def)),data=out.group3)+
  geom_point(aes(x=t,y=rel,shape="0.1",colour=Co_r),data=out.group1.thin)+
#  geom_point(aes(x=t,y=rel,shape="0.2",colour=Co_r),data=out.group2.thin)+
#  geom_point(aes(x=t,y=rel,shape="0.05",colour=Co_r),data=out.group3.thin)+
  #facet_grid(.~alpha)+
  geom_vline(xintercept=3,linetype=2,size=0.2)+
  geom_vline(xintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=16,linetype=2,size=0.2)+
  geom_text(aes(x=3, label='T[ps]', y=0.1), colour="blue", angle=0,parse=T)+
  geom_text(aes(x=10, label="T[h_s]", y=0.1), colour="red", angle=0,parse=T)+
  geom_text(aes(x=16, label="T[hs]+T[hd]", y=0.1),hjust=0.2, colour="red", angle=0,parse=T)+
  ylab("Relative difference in LAI (-)")+
  xlab("time (days)")+
  scale_y_continuous(limits=c(-0.05,0.1))+
  #scale_colour_manual("Co_r",values=col1)+
  scale_shape_manual("fPrim",values=points)+
  theme(legend.position="right",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  guides(colour= guide_legend(title.position="top",override.aes = list(shape = NA)))
dev.off()


##############################################
########### sensitivity analysis  ############
##############################################


#### Dynamics LAI
# set parameters
fPrim = c(0)
C_def = c(7.2)
D = c(8)
Ths = c(10)
The = c(180)
Tps = c(-7)
C_prim =  0.1#c(0.1,0.2)
Co_r = c(1*1e-08)
Thd = 6
Tpd = 100
alpha = c(1)
feeding = c("lognorm")
responsepat = c("Ea")

# make parameter matrix
a = expand.grid(C_def=C_def,D=D,Th_s= Ths,Th_e =The,Tp_s = Tps,Co_r=Co_r,fPrim=fPrim,C_prim=C_prim,Th_d = Thd,Tp_d=Tpd,alpha=alpha,responsepat=responsepat,feeding=feeding)
a$Tp_ss = a$Tp_s
a$Tp_dd = a$Tp_d
a$Tp_s = pmax(1,a$Th_s + a$Tp_s )
a$Tp_d = pmin(a$Th_s, a$Tp_d + a$Tp_s) - a$Tp_s
a = unique(a[,-14])


# add to full parameter dataframe
comb = read.csv("M:/My Documents/Plantcommunication/20160412_Input_defense_evo_Priming_onoff.csv")
pars = as.data.frame(matrix(NA,nrow=nrow(a),ncol=length(comb[1,])))

names(pars) = names(comb)
pars[] = comb[1,]

#a$C_prim = 1/a$C_def



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
pars$C_def_max = pars$C_def
pars$Ts = 180
pars$no = seq(1:nrow(pars))
pars$mut.val = 1
pars$Tp_d = a$Tp_d
pars$Th_d = 8
pars$Tco_d = 6
pars$shape = "linear"
pars$alpha = a$alpha
pars$feeding = a$feeding
pars$responsepat = a$responsepat
pars$Ft = 80
pars$feeding = "lognorm"

# parameters for which sensitivity analysis is done
sens.pars = c("Th_s","Tp_s","C_def","C_prim","Co_r","Tp_d","Th_d","Tco_d","alpha")

# set reference simulations
pars$sim = "ref"
pars$parname = "ref"
pars$C_prim = pars$C_prim_f
pars$fPrim = pars$fPrim_f
# add changed params for data.frame pars
for (i in 1:(length(sens.pars))){
  # add lower
  sens = pars[1,]
  sens$parname = sens.pars[i]
  sens$sim = "min"
  sens[,sens.pars[i]] = sens[,sens.pars[i]]*0.9
  pars = rbind(pars,sens)
  # add upper
  sens = pars[1,]
  sens$parname = sens.pars[i]
  sens$sim = "max"
  sens[,sens.pars[i]] = sens[,sens.pars[i]]*1.1
  pars = rbind(pars,sens)
}

pars$C_prim_f  = pars$C_prim
pars$C_prim_r  = pars$C_prim
pars$fPrim_f = pars$fPrim
pars$fPrim_r = pars$fPrim
pars$C_def_f = pars$C_def
pars$C_def_r = pars$C_def
pars$C_def_max = pars$C_def_f
pars = pars[,-c(47,48)]

results1 = find_endpoints_defense(pars,Ft=F)
results1 = data.frame(pars,results1)
results1$Tp_ss_1 = -(results1$Th_s - results1$Tp_s)
# store results
#write.csv(results1[[1]],paste(format(Sys.time(), "%Y%m%d"),"_output_Prim_GtOpt_Attack_NonAttack_linear8_delay.csv",sep=""))
write.csv(results1,paste(format(Sys.time(), "%Y%m%d"),"_Sims_Prim_Attack_NonAttack_sim_8_linear_delay6_sensanalysis.csv",sep=""))



sensitivity = function(x,results1){
  ref = results1[1,"seeds_tfm"]
  min = x[x$sim =="min","seeds_tfm"]  
  max = x[x$sim =="max","seeds_tfm"]  
  out = abs(c(((min-ref)/ref)*100,((max-ref)/ref)*100))
  #print(out)
}

sensitivity(results2[[1]],results1)

results2 = results1[-1,]
results2 = split(results2,results2$parname)
sens = t(sapply(results2,sensitivity,results1=results1))

apply(sens,1,mean)


### analysis Chris Frost

ir = read.csv("20170721_par.var.e.par.var.e.responsepat_Dvar_2.csv")

# remove points where optimal priming level is 0

# 
ir$rel. = (-ir$rel*ir$seeds_tf - -ir$ir*ir$seeds_tf)/(-ir$ir*ir$seeds_tf)

ir$Th_d = as.factor(ir$Th_d)
ir$D = as.factor(ir$D)
names(ir)[2] = "C_def"
ir.melt = melt(ir[,c(1,2,3,6,5,8)],id.var=1:4)

ir.melt = ir.melt[!(ir.melt$fPrim <0.0001 & ir.melt$variable=="rel"),]


ggplot(data=ir)+
  geom_line(aes(x=C_def,y=fPrim,colour=D,linetype=Th_d))

ir.1 = ir#[ir$D == 1e-08,]

]

pdf("inducible defense.pdf",width=6,height=4)
ggplot(data=ir.melt)+
  geom_line(aes(x=C_def,y=-value,colour=Th_d,linetype=variable))+
  facet_grid(.~ D)+
  theme(legend.position="right",
        axis.text.x=element_text(size=9), axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),axis.title.x=element_text(size=11),
        panel.background = element_rect(fill = "white", 
        colour = NA), panel.border = element_rect(fill = NA, 
        colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
        size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
        size = 0.5), strip.background = element_rect(fill = "grey80", 
        colour = "grey50"), legend.text=element_text(size=11),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9))+
  
  labs(linetype="Strategy",colour=expression(paste(italic(t)[hd])))+
  ylab("Relative difference in seed production (%)")+
  xlab(expression(paste(italic(C)[def])))
dev.off()
  
  #geom_line(aes(x=value,y=-ir,colour=D,linetype=Th_d))

ggplot(data=ir)+
  geom_line(aes(x=C_def,y=rel.,linetype=Th_d,colour=D))
  


