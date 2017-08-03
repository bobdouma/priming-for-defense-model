# Flowering and herbivory, priming and defense model 
# model of Peter Vermeulen (peter1.vermeulen@wur.nl) published in Vermeulen, New Phytoloist 2015 205:429-439
# rewritten to R by Bob Douma for the purpose of the priming paper (bob.douma@wur.nl). 




# version 2: includes the effect of herbivory on plant growth. Version tested and works fine now. 
#--- insect herbivore eats with constant rate a certain amount of leaf area per unit time
#--- the plant produces secondary compounds to decrease herbivory. This translates to higher construction costs. I assume that this does not affect leaf respiration, nor photosynthesis rates
#--- herbivory starts and ends at a given point in time

# version 3: as version 2. In addition: 
#--- includes a moment in time when emitter plant emits volatiles: Tp_s
#--- from that moment onwards a plant moves in a primed state (if Prim.r/Prim.f=T)
#--- The metabolic costs of priming are a fraction C_prim of the costs of defense
#----Reduction of consumption rate are a fraction of consumption rate when defended  when plant is primed; fPrim. 
#----Maintenance costs scale linearly with C_def (that if leave is defended; maintenance costs scale proportionally with C_def (as construction costs))

# this version as previous: except:
#--- the plant produces secondary compounds to decrease herbivory. This translates to higher maintenance costs. I assume that construction costs are incorporated in the maintenance costs
#--- the relationship between C_def and CO_r is made linear

# version 2: as previous version:  
#--- defense function has different shape. The reduction in defense can have a lag time if needed
#--- priming only takes place between the release of the cue and the actual attack. Afterwards the plant pays defense costs, but with no lag time

# version 3: as previous version:
#--- feeding rate increases in X days form zero to maximum level. Parameter Tco_d added to indicate delay in feeding rate

# verion 4: as previous version:
#--- priming is the first part of defense; during defense only the part of defense needs to be build up that was not build up during priming 

# version 5: as previous versions
#--- feeding rate increases in 12 days form zero to maximum level. Parameter Tco_d added to indicate delay in feeding rate
#----competition can be asymmetric depending on the symmetry factor alpha


library(utils)
library(compiler)
library(ggplot2)
library(foreach)  
library(doParallel)

# remove objects from workspace
#rm(list=ls(all=TRUE))
# set working directory
#setwd("D:\\BobDouma\\Veni\\Game_theory_plantcommunication\\simulations\\sim_8_linear\\")


##### Initial parameters #######
## constants etc

#time
Ts= 180.0  #season length in days
h= 24.0  #hours in a day
I_day= h #required intervals per day
D_h= h/2  #day light hours

#photosynthesis
k= 0.8  #light extinction coefficient 
Lue_st= 0.020 #mol mol-1; Hikosaka et. 1999 Oecologia; taken as average high Lue; note in the paper, Lue changes with plant mass
Io= 1000.0  #umol quanta s-1

#allocation parameters
F_la= 0.33  #fraction of mass to leaf mass, assuming equal sopport of rootsd and stem
S_m= 0.01  #seed weight in g
SLA= 0.03  # m2 g-1 Pronk 2004

#respiration constants
C_built_resp= 0.213  #gC g-1C into leaf  #respiration rate in C to built 1g of C into new leaves
Leaf_c= 0.45  #gC g-1 leaf; so gram C in 1 gram of leaf Pronk 2004

##calculation of constants
S_day= D_h*60*60 #Seconds per day of light

# respiration and leaf cost functions
# building costs and building respiration in terms umol C + costs of other tissues to support leaves (1/F_la)
C_newleaf= (1/F_la)*((Leaf_c*1.0+ Leaf_c*C_built_resp)/12.011)*1000000.0  #umol C g-1 leaf; resp gC per gram g leaf built / molar weight of C*conversion from mol to umol
C_def = 8.71747 # umol m-2 s-1  costs invested in defense; calculated from van Hulten et al 2006 PNAS (calculation based on reduction in RGR)
# not a parameter but calculated from van hulten 2006 PNAS
Cost_prim = 3.359094 # umol m-2 s-1  costs of priming (calculated from van Hulten et al 2006 PNAS and solving with Basic LAI model)
C_prim = Cost_prim/C_def # fraction of priming costs compared to defense

# Consumption rate of herbivores
Co_r = 0.0000002 #consumption rate of one insect herbivore  (m2/s)
Co_r = 8e-09 # m2/s average consumption rate of one P.rapae over its lifetime
fPrim = 0.06 # factor that modulates the reduction in consumption rate (Co_r.r). If primed, reduction is stronger compared to non-primed

F_nightresp= h/D_h  #accounts for night time respiration as a ratio of daylength/ light length
C_l= 1.166 #umol m-2 s-1 #calculated from Pronks maintenance constant in gC/gC mass/ day, seconds per day of 25200, SLA, and C content leaf= 0.45
# (1/F_la) costs of other organs are included in the Pnet.
C=C_l*(1/F_la)  # as a percentage of Pmax* a correction factor that accounts for night time respiration; set h to D_h when assuming no night time respiration


#### helper functions

# inverse logit
inv.logit = function(x){
  exp(x)/(exp(x)+1)
}

# function to extract a column from a dataframe with name tat is stored in mod
fun.name = function(x,mod){
  a = x[,eval(mod)]
  return(a)
}


#function to calculate start LAI from seed mass, density, the fraction of mass to leaf area, and SLA
Start_LAI = function(S_m, d, F_la, SLA){  
  LAI= S_m*d*F_la*SLA  #seed mass (g seed-1)*density of plant (# seeds / m2 ground surface)*fraction to leaf area* leaf area per unit mass (m2 g-1)
  return(LAI)
}


# function to calculate gross photosynthesis
Pgross = function(Io, k, LAI_t, Lue){
  Pg= Lue*Io*(1-exp(-k*LAI_t))   #Monsi and Saeki light intercepted by canopy*light use efficiency
    return(Pg) ##equal to Io-Io*e-k*lai; outcome in umol Co2 s-1 m-2 area soil surface
}

# function to calculate Pnet for focal plant with leaf area_f in competition with resident plants with Leaf area_r
Pnet = function(LAI_f, LAI_t, Pgr,D, D_r, C, C_def,C_def_max, C_prim,t,Th_s.i,Th_e.i,Tp_s.i,Tp_e.i,alpha,fPrim,Herbivory, Defense, Priming){ 
  if (t < (Th_s.i) | t > (Th_e.i)) {Herbivory = F; Defense= F} # assume that defense is off when herbivory is absent. t is in numbers (1:Ts*I_day) 
  if (Priming & (t > Tp_s.i)){C_def = min(C_def_max,C_def*(1+Priming*fPrim))} # Higher investment in defense when primed if fPrim > 1
  if (Priming & (t < Tp_s.i | t > Tp_e.i)) {Priming=F} # assumes that priming equals defense. Priming stops when herbivory stops
  if (Pgr == 0){Pnet_f = 0} else { # if gross photosynthesis is zero Pnet is 0, otherwise:
    if (D_r == 0){fr = 1} else if (D ==0) {fr =0} else {
      fr = (D*((LAI_f/D)/(LAI_t/(D_r + D)))^alpha)/((D*((LAI_f/D)/(LAI_t/(D_r + D)))^alpha)+(D_r*(((LAI_t-LAI_f)/(D_r))/(LAI_t/(D_r + D)))^alpha))# competition coefficient
    }
    Pnet_f= fr*Pgr- (LAI_f* C) - (C_def* Defense*LAI_f) - (C_prim*C_def*Priming*LAI_f) # assumes  that maintenance costs (also with respect to priming and defense) scale proportionally with LAI. The C_def relative to C_def_max determines 
  }  #fraction capture by focal plant* Pgross of vegetation - respiration costs; C= cost per m2 -2 leaf
  return(Pnet_f) #outcome is Pnet of focal plant (in umol co2 s-1 m-2 soil surface)
}


# function to calculate  for t =FT what the total seed output over the season would be.
Seeds_season_simple = function(Pn, Ts, t){ 
  #so Pn is net photosynthesis of the focal plant, Ts is length of season, and t is the time of interest
  Seeds= Pn*(Ts-t) #not that Pnet is now in umol CO2 -s, and so Ts and T should also be in seconds
  return(Seeds)
}

## This calculates the feeding rate of the insect herbivore based on the level of defense and the delay in activation defenses
feedingrate = function(C_def,C_def_max,C_prim,Co_r,dt,t,Th_s.i,Th_e.i,Tp_s.i,Tp_e.i,Th_d.i,fPrim,shape,responsepat,feeding,Herbivory, Defense, Priming){
  if (t < (Th_s.i) | t > (Th_e.i)) {Herbivory = F; Defense= F} # assume that defense is off when herbivory is absent. t is in numbers (1:Ts*I_day) 
  if (responsepat == "Fa"){Fa=T;Ea=F; St=F} else if (responsepat == "Ea"){Fa=F;Ea=T;St=F} else if (responsepat == "St"){Fa=F;Ea=F;St=T} else {stop("responsepat not specified")}
  
  if (Priming){
    Th_d.i = Th_d.i+(Th_d.i-Th_s.i)*(min((C_def_max/C_def),fPrim)*Priming) # if C_def is fPrim higher, also the time to build up defense is higher 
    # time that primed plant is ahead of defended plant in terms of building up defense
    t.advance.b = min(C_prim*(Th_d.i-Th_s.i), (Th_s.i-Tp_s.i))
    #t.advance.b = min(0.5*(Th_d.i-Th_s.i), (Th_s.i-Tp_s.i))
    # time since plant was in primed state
    t.advance.d = max(0,C_prim*(Th_d.i-Th_s.i)-(Th_s.i-Tp_e.i))
    #t.advance.d = max(0,0.5*(Th_d.i-Th_s.i)-(Th_s.i-Tp_e.i))
    # determines the effective time priming is ahead of non-primed plant
    t.advance = min(t.advance.b,t.advance.d)
  } else {t.advance = 0; Fa=F;Ea=F;St=F} # # non-primed has no time advantage in buildup of defenses
  
  C_def = min(C_def_max,C_def*(1+fPrim*St)) # when priming responsepatter is "stronger" it builds up a higher defense level
  # calculate actual defense level depending on the responsepattern of a  primed plant and the t.advance. If priming is off there is a delay in the build up of defense determined by Th_d.i - Thsi
  rC_def = min((C_def),((C_def)/(Th_d.i-t.advance*Fa-Th_s.i))*((t+t.advance*Ea)-Th_s.i),na.rm=T) #  level of defense corrected start up time of defense (+dt to prevent dividing by zero)
  # calculates effective reduction in herbivory based on assumtpion that
  # feeding rate reduces linearly with defense level
  if (shape=="linear"){
    Co_r.r =  min(1,(1/C_def_max)*rC_def) #  linear reduction in herbivory by increase in defense - varies between 0 and 1  by scaling to C_def_max
  } else if (shape =="sigmoid"){
    Co_r.r = inv.logit(-4.59+((4.59--4.59)/C_def_max)*rC_def) #  sigmoidal reduction in herbivory by increase in defense (at defense level max, reduction =0.99, at defense =0, reduction =0.01)
  } else if (shape == "delay"){print("not implemented anymore")}
  
    
  rCo_r = Co_r * (1-(Co_r.r*Defense)) # calculates the reduced Co_r based on the potential feeding rate and the level of defense
    
  return(rCo_r)  
}




# function to calculate the net change in LAI, to calcualte total LAI at next time step
LAI_delta = function(D,Pn, LAI, C_newleaf,SLA,C_def,C_def_max,C_prim,Co_r,dt,t,Th_s.i,Th_e.i,Tp_s.i,Tp_e.i,Th_d.i,Tco_d.i,fPrim,responsepat,shape,feeding,Herbivory, Defense, Priming){
  if (t < (Th_s.i) | t > (Th_e.i)) {Herbivory = F; Defense= F} # assume that defense is off when herbivory is absent. t is in numbers (1:Ts*I_day) 
  # the leaf area consumed by insect herbivore (either for focal or resident plant)
  ####### three ways to specify consmption rate

  if (Herbivory){
    if (feeding == "constant"){
      # 1 linear increase until maximum then constant
      Co_r = D * max(min(Co_r/Tco_d.i*(t-Th_s.i),Co_r,na.rm=T),0)   # D plant density   
    } else if (feeding =="lognorm"){
      # 2 following community dynamics
      Co_r = D * Co_r*exp(-(log((t-Th_s.i)/24)-3.39)^2/(2*0.6167^2))  ## mean: 7.67*exp(-(log((t-Th_s.i)/24)-3.39)^2/(2*0.6167^2)); # D plant density   
    } else if (feeding =="prop") {
      # 3 proportion of LAI
      Co_r = Co_r * LAI
    } else {stop("feeding not specified")}
    # reduction in feeding rate
    rCo_r = feedingrate(C_def=C_def,C_def_max=C_def_max,C_prim=C_prim,Co_r=Co_r,dt=dt,t=t,Th_s.i=Th_s.i,Th_e.i=Th_e.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Th_d.i=Th_d.i,fPrim=fPrim,shape=shape,responsepat=responsepat,feeding=feeding,Herbivory=Herbivory, Defense=Defense, Priming=Priming)
    d_LAI_h = (Herbivory * rCo_r) #m2 # Co_r = consumption rate   (m2/s), Co_r.r = reduction in consumption rate (-); D = number of plants per m2
  } else {Co_r =0; d_LAI_h =0; rCo_r = Co_r}
  #  if (t > 970 & t< 980){
  #  x = data.frame("ts"=t,"t"=t/43200,"rCo_r"=rCo_r[1],"Co_r"=Co_r,"C_def"=C_def,"rC_def"=rCo_r[2],"Co_r.r"=rCo_r[3],"Thsi"=Th_s.i,"Thdi"=Th_d.i,"d_LAI_h"=d_LAI_h,"Her"=Herbivory,"Def"=Defense,"Priming"=Priming)
  #  write.table(x,"output.csv",col.names=FALSE,append=TRUE,sep="\t")
  # }
  # the construction costs of the leaf (defended or not)

  PtoLa = C_newleaf  /SLA # building costs and building respiration in terms unit area m-2. Default costs C_newleaf, and additional costs, which is a fraction of the default costs 
  d_LAI= (Pn/PtoLa-d_LAI_h)*dt #Pn= net photosynthesis (umo)l, PtoLA is a conversion factor for net P to umol C m-2 leaf area m-2 soil surface; we still need time; so here dt is in seconds (see dt_s in euler); 
  #  print(c(C_prim,C_def,Defense,Co_r.r,d_LAI_h,PtoLa,D,Herbivory,Co_r))
  return(d_LAI)
}


# function that calculates the benefit of a mutant at LAI t to delay flowering by dt
Seeds_m = function(LAI_f, LAI_r, dLAI,D,D_r, Ts, t, dt,k,Lue,C,Io, C_def, C_def_max, C_prim,fPrim,alpha,Th_s.i,Th_e.i,Tp_s.i,Tp_e.i,Herbivory, Defense, Priming){  
  LAI_mu= max(0,LAI_f+ dLAI)  #LAI of mutant at time t+dt
  LAI_t= LAI_mu+ LAI_r  #LAI total at time t+dt of the whole vegetation
  Pgr= Pgross(Io=Io, k=k, LAI_t=LAI_t, Lue=Lue)  #consequent Pgross at time t+dt of whole vegetation
  Pnet_f= Pnet(LAI_f=LAI_mu, LAI_t=LAI_t, Pgr=Pgr, C=C,D = D, D_r = D_r, C_def=C_def,C_def_max=C_def_max,C_prim=C_prim,fPrim=fPrim, alpha=alpha,t=t,Th_s.i=Th_s.i, Th_e.i=Th_e.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i, Herbivory=Herbivory, Defense=Defense, Priming=Priming)
  
  # if Pn < 0 during seed production (because of high defense costs); seed production is not possible; Pn is set to zero
  if (Pnet_f < 0){Pnet_f =0}
  
  Seed_mu= Seeds_season_simple(Pn=Pnet_f, Ts= Ts, t = (t))  #seed production if flowering time of mutant= t+dt, while FT of resident remains at t
  # from Seed_mu removed dt from t = t+dt
  #delta_seed_m= (Seed_tm1-Seeds_f)/dt #note we compare fitness against the focal plant of equal frequency as the mutant
  #so, we assume that the resident population is so large that the mutant does not affect its fitness; also note it calcuates the slope (dividing by dt)    
  return(Seed_mu)
}

# calculate the seed production over the season
Seeds_season = function(no,D_r,Pn_f,Pn_r, Ts,dt,t,j,n,LAI_f0,LAI_r0,dLai_tf,dLai_tr,Th_e.i, Th_s.i,Th_d.i,Tp_s.i,Tp_e.i,Tco_d.i,
                        C_def_f,C_def_r,C_def_max,C_prim_f,C_prim_r,C,C_newleaf,Co_r, fPrim_f,fPrim_r,alpha,
                        SLA,k,Lue,S_day,Io,H.r,H.f,Def.r,Def.f,Prim.r,Prim.f,shape,responsepat,feeding){   # t comes in as a vector with seconds
  #Ts is in seconds
  LAI_f = LAI_f0 # start LAI values 
  LAI_r = LAI_r0  # start LAI values 
  #  LAI_mu= max(0,LAI_f+ dLAI) #LAI of mutant at time t+dt
  
  # if Pn < 0 during seed production (because of high defense costs); seed production is not possible; Pn is set to zero
  if (Pn_f < 0){Pn_f =0}
  if (Pn_r < 0){Pn_r =0}
  # store output
  n1 = 180*24
  LAI_f_o= rep(0,n1)
  LAI_r_o= rep(0,n1)
  dLAI_tf_o = rep(0,n1)
  dLAI_tr_o = rep(0,n1)
  LAI_mu_f_o = rep(0,n1)
  LAI_mu_r_o = rep(0,n1)
  Pn_f_o= rep(0,n1)
  Pn_r_o= rep(0,n1)
  Pn_mu_f_o= rep(0,n1)
  Pn_mu_r_o= rep(0,n1)
  seeds_tf_o = rep(0,n1)
  seeds_tr_o = rep(0,n1)
  seeds_mu_o = rep(0,n1)
  
  # seed production for the day at which plant starts flowering
  seeds_tf = Pn_f*dt # seed production of focal for t =j
  seeds_tr = Pn_r*dt # seed production of resident t =j
  seeds_mu = 0 # seed production of mutant that delays flowering with one time step
  seeds_tf_c = seeds_tf
  seeds_tr_c = seeds_tr
  seeds_mu_c = seeds_mu
  # three situations are distinguished:
  # 1: when flowering starts before herbivory is started
  # 2: when flowering starts when herbivores are present
  # 3: when flowering starts when herbivores are not present (and will not come)
  # situation :
  
  if ((t[j+1]) < t[Th_s.i] & t[j+1] < Ts){
    #    print("first loop")
    # before herbivory
    seeds_tf = seeds_tf + Seeds_season_simple(Pn = Pn_f,Ts=(t[Th_s.i]),t=(t[j+1]))
    seeds_tr = seeds_tr + Seeds_season_simple(Pn = Pn_r,Ts=(t[Th_s.i]),t=(t[j+1]))
    seeds_mu = seeds_mu + Seeds_m(LAI_f = LAI_f,LAI_r = LAI_r,dLAI = dLai_tf,D=1,D_r = D_r,Ts=(t[Th_s.i]),t=(t[j+1]),dt=dt,k=k,Lue=Lue,C=C,Io=Io,alpha = alpha,
                                  C_def_f =C_def_f,C_def_max=C_def_max,C_prim=C_prim_f,fPrim=fPrim_f,Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f,Defense=Def.f,Priming=Prim.f)
    
    t_h = t[(t >= t[Th_s.i] & t < t[Th_e.i])] # vector with t seconds in which herbivore is active
    #    print(length(t_h))
    #    print(Ts)
    LAI_mu = dLai_tf + LAI_f
    # during herbivory
    for (i in 1:length(t_h)){ # seed production for the remaining season 
      
      dLai_tf= LAI_delta(D=1,Pn=0, LAI = LAI_f, C_newleaf = C_newleaf,SLA = SLA,C_def_f = C_def_f, C_def_max=C_def_max, C_prim=C_prim_f,Co_r = Co_r,fPrim = fPrim_f,dt = dt,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Th_d.i=t[Th_d.i],Tco_d.i=t[Tco_d.i], Herbivory=H.f,Defense=Def.f,Priming=Prim.f,shape=shape,responsepat=responsepat,feeding=feeding)  # if Pnet goes to LAI, this is the increase in LAI over dt; focal plant
      dLai_tr= LAI_delta(D=D_r,Pn=0, LAI = LAI_r, C_newleaf = C_newleaf,SLA = SLA,C_def_r = C_def_r, C_def_max=C_def_max, C_prim=C_prim_r,Co_r = Co_r,fPrim = fPrim_r,dt = dt,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Th_d.i=t[Th_d.i],Tco_d.i=t[Tco_d.i], Herbivory=H.r,Defense=Def.r,Priming=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)  #if Pnet goes to LAI, this is the increase in LAI over dt; resident plants
      LAI_f = max(0,(LAI_f + dLai_tf))  # calculate new leaf area of focal
      LAI_r = max(0,(LAI_r + dLai_tr))  # calculate new leaf area of resident
      LAI_t= LAI_f + LAI_r # calcalute total leaf area
      Pgr= Pgross(Io, k, LAI_t, Lue)  #Gross photosynthesis of the whole vegetation; with LAI_total= LAI_focal+ LAI resident
      Pn_f= Pnet(LAI_f = LAI_f, LAI_t = LAI_t, Pgr = Pgr, D = 1, D_r = D_r , C= C,C_def =C_def_f,C_def_max=C_def_max,C_prim=C_prim_f,fPrim=fPrim_f,alpha = alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f,Defense=Def.f,Priming=Prim.f)  # net photosynthesis as a share of Pgross (ratio of LAI? LAI_t) - the cost of the focal plant
      Pn_r= Pnet(LAI_f = LAI_r, LAI_t = LAI_t, Pgr = Pgr, D = D_r, D_r = 1, C= C,C_def =C_def_r,C_def_max=C_def_max,C_prim=C_prim_r,fPrim=fPrim_r,alpha = alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.r,Defense=Def.r,Priming=Prim.r)  # net photosynthesis as a share of Pgross (ratio of LAI/LAI_t) - the cost of the resident plants
      # if Pn < 0 during seed production (because of high defense costs); seed production is not possible; Pn is set to zero
      if (Pn_f < 0){Pn_f =0}
      if (Pn_r < 0){Pn_r =0}
      if(Pn_f < 0){print("loop 1")}
      seeds_tf = seeds_tf + Pn_f*dt # seed production of focal in this time step
      seeds_tr = seeds_tr + Pn_r*dt # seed production of resident in this time step
      # seed production of mutant that flowers dt later (i.e. has to be compared to Seeds_tf before the for loop)
      LAI_mu = max(0,dLai_tf + LAI_mu) 
      LAI_t = LAI_mu  + LAI_r  ### adapt!!! #LAI total at time t+dt of the whole vegetation
      Pgr= Pgross(Io=Io, k=k, LAI_t=LAI_t, Lue=Lue)  #Gross photosynthesis of the whole vegetation; with LAI_total= LAI_focal+ LAI resident
      Pn_mu= Pnet(LAI_f = LAI_mu, LAI_t=LAI_t, Pgr=Pgr, D = 1, D_r = D_r, C=C,C_def =C_def_f,C_def_max=C_def_max,C_prim=C_prim_f,fPrim=fPrim_f, alpha = alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f,Defense=Def.f,Priming=Prim.f) ### adapt # net photosynthesis as a share of Pgross (ratio of LAI? LAI_t) - the cost of the focal plant
      # if Pnet of mutant is < 0 no seed production is possible
      if (Pn_mu < 0){Pn_mu =0}
      seeds_mu = seeds_mu + Pn_mu*dt
      # seed production of focal and resident later in season. First calculate new leaf area due to possible herbivory
      
    }
    # after herbivory
    seeds_tf = seeds_tf + Seeds_season_simple(Pn = Pn_f,Ts=Ts,t=t[Th_e.i])
    seeds_tr = seeds_tr + Seeds_season_simple(Pn = Pn_r,Ts=Ts,t=t[Th_e.i])
    seeds_mu = seeds_mu + Seeds_m(LAI_f = LAI_mu,LAI_r = LAI_r,dLAI = 0,D = 1, D_r = D_r, Ts=Ts,t=t[Th_e.i],dt=dt,k=k,Lue=Lue,C=C,Io=Io,
                                  alpha = alpha, C_def =C_def_f,C_def_max=C_def_max, C_prim=C_prim_f,fPrim=fPrim_f,Th_s.i= t[Th_s.i],Th_e.i = t[Th_e.i],Tp_s.i = t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f, Defense=Def.f, Priming=Prim.f)
    #    print(c("seedsmu",seeds_mu))
    
  }   else if (((t[j+1]) <= t[Th_e.i]) & ((t[j+1]) >= t[Th_s.i])){
    #    print("second loop")
    t_h = t[(t >= (t[j+1]) & t <= t[Th_e.i])] # vector with t seconds in which herbivore is active
    LAI_mu = dLai_tf + LAI_f
    for (i in 1:length(t_h)){ # seed production for the remaining season 
      
      # seed production of focal and resident later in season. First calculate new leaf area due to possible herbivory
      dLai_tf= LAI_delta(D=1,Pn=0, LAI = LAI_f, C_newleaf = C_newleaf,SLA = SLA,C_def = C_def_f,C_def_max=C_def_max,C_prim=C_prim_f,Co_r = Co_r,fPrim = fPrim_f,dt = dt,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Th_d.i=t[Th_d.i],Tco_d.i=t[Tco_d.i], Herbivory=H.f,Defense=Def.f,Priming=Prim.f,shape=shape,responsepat=responsepat,feeding=feeding)  # if Pnet goes to LAI, this is the increase in LAI over dt; focal plant
      dLai_tr= LAI_delta(D=D_r,Pn=0, LAI = LAI_r, C_newleaf = C_newleaf,SLA = SLA,C_def = C_def_r,C_def_max=C_def_max,C_prim=C_prim_r,Co_r = Co_r,fPrim = fPrim_r,dt = dt,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Th_d.i=t[Th_d.i],Tco_d.i=t[Tco_d.i], Herbivory=H.r,Defense=Def.r,Priming=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)  #if Pnet goes to LAI, this is the increase in LAI over dt; resident plants
      LAI_f = max(0,(LAI_f + dLai_tf))  # calculate new leaf area of focal
      LAI_r = max(0,(LAI_r + dLai_tr))  # calculate new leaf area of resident
      LAI_t= LAI_f + LAI_r # calcalute total leaf area
      Pgr= Pgross(Io, k, LAI_t, Lue)  #Gross photosynthesis of the whole vegetation; with LAI_total= LAI_focal+ LAI resident
      Pn_f= Pnet(LAI_f = LAI_f, LAI_t = LAI_t, Pgr= Pgr,D=1, D_r = D_r, C=C,C_def =C_def_f,C_def_max=C_def_max,C_prim=C_prim_f,fPrim=fPrim_f,alpha =alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f,Defense=Def.f,Priming=Prim.f)  # net photosynthesis as a share of Pgross (ratio of LAI? LAI_t) - the cost of the focal plant
      Pn_r= Pnet(LAI_f = LAI_r, LAI_t = LAI_t, Pgr = Pgr, D = D_r, D_r = 1, C= C,C_def =C_def_r,C_def_max=C_def_max,C_prim=C_prim_r,fPrim=fPrim_r,alpha =alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.r,Defense=Def.r,Priming=Prim.r )  # net photosynthesis as a share of Pgross (ratio of LAI/LAI_t) - the cost of the resident plants
      # if Pn < 0 during seed production (because of high defense costs); seed production is not possible; Pn is set to zero
      if (Pn_f < 0){Pn_f =0}
      if (Pn_r < 0){Pn_r =0}
      if(Pn_f < 0){print("loop 2")}
      seeds_tf = seeds_tf + Pn_f*dt # seed production of focal in this time step
      seeds_tr = seeds_tr + Pn_r*dt # seed production of resident in this time step
      
      LAI_mu = max(0,dLai_tf + LAI_mu) 
      # seed production of mutant that flowers dt later (i.e. has to be compared to Seeds_tf before the for loop)
      LAI_t = LAI_mu  + LAI_r  ### adapt!!! #LAI total at time t+dt of the whole vegetation
      Pgr= Pgross(Io=Io, k=k, LAI_t=LAI_t, Lue=Lue)  #Gross photosynthesis of the whole vegetation; with LAI_total= LAI_focal+ LAI resident
      Pn_mu= Pnet(LAI_f = LAI_mu, LAI_t=LAI_t, Pgr=Pgr,D=1,D_r = D_r, C=C,C_def =C_def_f, C_def_max=C_def_max,C_prim=C_prim_f,fPrim=fPrim_f,alpha =alpha,t=t_h[i],Th_s.i=t[Th_s.i],Th_e.i=t[Th_e.i],Tp_s.i=t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f,Defense=Def.f,Priming=Prim.f) ### adapt # net photosynthesis as a share of Pgross (ratio of LAI? LAI_t) - the cost of the focal plant
      # if Pnet of mutant is < 0 no seed production is possible
      if (Pn_mu < 0){Pn_mu =0}
      seeds_mu = seeds_mu + Pn_mu*dt
      # store output
      LAI_f_o[i] = LAI_f
      LAI_r_o[i] = LAI_r
      dLAI_tf_o[i] = dLai_tf
      dLAI_tr_o[i] = dLai_tr
      LAI_mu_f_o[i] = LAI_mu
      Pn_f_o[i]= Pn_f
      Pn_r_o[i]= Pn_r
      Pn_mu_f_o[i] = Pn_mu
      seeds_tf_o[i] = seeds_tf
      seeds_tr_o[i] = seeds_tr
      seeds_mu_o[i] = seeds_mu
    }
    # after herbivory
    seeds_tf = seeds_tf + Seeds_season_simple(Pn = Pn_f,Ts=Ts,t=t[Th_e.i+1])
    seeds_tr = seeds_tr + Seeds_season_simple(Pn = Pn_r,Ts=Ts,t=t[Th_e.i+1])
    seeds_mu = seeds_mu + Seeds_m(LAI_f = LAI_mu,LAI_r = LAI_r,dLAI = 0,D=1,D_r=D_r,Ts=Ts,t=t[Th_e.i+1],dt=dt,k=k,Lue=Lue,C=C,Io=Io,alpha =alpha,
                                  C_def =C_def_f,C_def_max=C_def_max, C_prim=C_prim_f,fPrim=fPrim_f,Th_s.i= t[Th_s.i],Th_e.i = t[Th_e.i],Tp_s.i = t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f, Defense=Def.f, Priming=Prim.f)
    #    }
  } else {
    #      print("third loop")
    # after herbivory
    if(Pn_f < 0){print("loop 3")}
    seeds_tf = seeds_tf + Seeds_season_simple(Pn=Pn_f,Ts=Ts,t=(t[j+1]))
    seeds_tr = seeds_tr + Seeds_season_simple(Pn=Pn_r,Ts=Ts,t=(t[j+1]))
    seeds_mu = seeds_mu + Seeds_m(LAI_f=LAI_f,LAI_r = LAI_r,dLAI= dLai_tf,D=1,D_r=D_r,Ts=Ts,t=(t[j+1]),dt=dt,k=k,Lue=Lue,C=C,Io=Io,alpha =alpha,
                                  C_def =C_def_f, C_def_max=C_def_max, C_prim=C_prim_f,fPrim=fPrim_f,Th_s.i= t[Th_s.i],Th_e.i = t[Th_e.i],Tp_s.i = t[Tp_s.i],Tp_e.i=t[Tp_e.i],Herbivory=H.f, Defense=Def.f, Priming=Prim.f)
    LAI_mu = max(0,dLai_tf + LAI_f)
  }
  #  seeds_tf_c = seeds_tf_c + Seeds_season_simple(Pn=Pn_f_c,Ts=Ts,t=(t[j+1]))
  #  seeds_tr_c = seeds_tr_c + Seeds_season_simple(Pn=Pn_r_c,Ts=Ts,t=(t[j+1]))
  #  seeds_mu_c = seeds_mu_c + Seeds_m(LAI_f=LAI_f_c,LAI_r = LAI_r_c,dLAI= dLai_tf_c,Ts=Ts,t=(t[j+1]),dt=dt,k=k,Lue=Lue,C=C,Io=Io,
  #                                    C_def =C_def, C_prim=C_prim_f,Th_s.i= t[Th_s.i],Th_e.i = t[Th_e.i],Tp_s.i = t[Tp_s.i],Herbivory=H.f, Defense=Def.f, Priming=Prim.f)
  
  delta_seed_m= (seeds_mu-seeds_tf)/dt
  out = c(seeds_tf, seeds_tr,seeds_mu, delta_seed_m, LAI_mu,LAI_f,LAI_r,dLai_tf,Pn_f,dt,t[j])
  names(out) = c("seeds_tf","seeds_tr","seeds_mu","delta_seed_m","LAI_mu","LAI_f","LAI_r","dLAI","Pn_f","dt","t" )
  out.sim = data.frame(LAI_f_o,LAI_r_o,dLAI_tf_o,dLAI_tr_o,Pn_f_o,Pn_r_o,seeds_tf_o,seeds_tr_o,seeds_mu_o)
  write.csv(out.sim,paste("Output_defense_",no,"_Defense_linear_8_delay_seeds_",format(Sys.time(), "%Y%m%d"),".csv",sep=""))
  
  return(out)
}



#function were we need to return several delta functions: delta lai of course, but also total seeds, and mutant values; so 3 outputs. 
outcomes_t = function(D_r, LAI_f, LAI_r, Io, k, Lue, C, C_newleaf,C_def_f,C_def_r,C_def_max,C_prim_f,C_prim_r,SLA, Co_r,fPrim_f,fPrim_r, dt, Ts,alpha,
                      t,Th_s.i,Th_e.i,Tp_s.i,Tp_e.i,Th_d.i,Tco_d.i,H.r,H.f,Def.r,Def.f,Prim.r,Prim.f,mod,shape,responsepat,feeding){
  #then in the Euler metod, first state  the three outcomes, and use  only dLai to iterate
  LAI_t= LAI_f+ LAI_r
  Pgr= Pgross(Io, k, LAI_t, Lue)  #Gross photosynthesis of the whole vegetation; with LAI_total= LAI_focal+ LAI resident
  Pn_f= Pnet(LAI_f = LAI_f, LAI_t = LAI_t, Pgr= Pgr, C=C, D=1, D_r = D_r ,C_def= C_def_f, C_def_max = C_def_max,C_prim=C_prim_f, fPrim=fPrim_f,alpha =alpha,t=t,Th_s.i=Th_s.i,Th_e.i=Th_e.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Herbivory=H.f, Defense=Def.f, Priming=Prim.f)  # net photosynthesis as a share of Pgross (ratio of LAI? LAI_t) - the cost of the focal plant
  Pn_r= Pnet(LAI_f = LAI_r, LAI_t = LAI_t, Pgr= Pgr, C=C, D = D_r, D_r = 1, C_def= C_def_r, C_def_max = C_def_max, C_prim=C_prim_r, fPrim=fPrim_r,alpha =alpha,t=t,Th_s.i=Th_s.i,Th_e.i=Th_e.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Herbivory=H.r, Defense=Def.r, Priming=Prim.r)  # net photosynthesis as a share of Pgross (ratio of LAI/LAI_t) - the cost of the resident plants
  dLai_tf= LAI_delta(D=1,Pn=Pn_f,LAI = LAI_f, C_newleaf = C_newleaf,SLA = SLA,C_def = C_def_f, C_def_max = C_def_max,C_prim = C_prim_f,Co_r = Co_r, fPrim=fPrim_f, dt = dt,t=t,Th_s.i=Th_s.i,Th_e.i=Th_e.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,shape=shape,responsepat=responsepat,feeding=feeding, Herbivory=H.f,Defense=Def.f,Priming=Prim.f)  # if Pnet goes to LAI, this is the increase in LAI over dt; focal plant
  dLai_tr= LAI_delta(D=D_r,Pn=Pn_r,LAI = LAI_r, C_newleaf = C_newleaf,SLA = SLA,C_def = C_def_r, C_def_max = C_def_max,C_prim = C_prim_r,Co_r = Co_r, fPrim=fPrim_r, dt = dt,t=t,Th_s.i=Th_s.i,Th_e.i=Th_e.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,shape=shape,responsepat=responsepat,feeding=feeding, Herbivory=H.r,Defense=Def.r,Priming=Prim.r)  #if Pnet goes to LAI, this is the increase in LAI over dt; resident plants
  Lai_tf= max(0,LAI_f+ dLai_tf) #calculate total leaf area for next time step (t+dt) for focal plant; LAI can't get below zero
  Lai_tr= max(0,LAI_r+dLai_tr)  #calcaulte total leaf area for next time step (t+dt) for focal plant; LAI can't get below zero
  return(c(LAI_f = Lai_tf, LAI_r = Lai_tr,P_gr = Pgr,Pn_f = Pn_f,Pn_r = Pn_r, dLAI_tf =t(dLai_tf),dLAI_tr =t(dLai_tr)))
}

# forward Euler scheme (Langtangen; 2009; page 615); adapted to produce multiple outputs so the seeds produced when Ft=t
Explicit_Euler_defense = function(i,D_r, S_m, F_la, SLA,Lai_f0, Lai_r0, Ts, I_day, S_day, Io,alpha, 
                                  k, Lue, C, C_newleaf ,C_def_f,C_def_r,C_def_max,C_prim_f, C_prim_r, Co_r,fPrim_f,fPrim_r,
                                  Ft,Th_s,Th_e,Tp_s,Tp_d,Th_d,Tco_d,Def.f,Def.r,H.f,H.r,Prim.f,Prim.r,mod,mut.val,shape,responsepat,feeding){ 
  
  # assign values to mutant (default eqaul to focal plant, except the parameter that is mentioned in mod)
  C_def_f_mu = C_def_f
  C_prim_f_mu = C_prim_f
  Co_r_mu = Co_r
  fPrim_f_mu = fPrim_f
  Def.f_mu = Def.f
  Prim.f_mu = Prim.f
  # modify value  of mutant that is entered in mod
  mod_mu = paste(mod,"_mu",sep="")
  
  # if Prim.f (T/F) then mut.val is added instead of multiplication
  if (eval(mod) == "Prim.f"){assign(mod_mu,get(mod_mu)+mut.val)} else {
    assign(mod_mu,get(mod_mu)*mut.val)}
  #print(c(mod,mod_mu,mut.val,get(mod_mu)))
  n= ceiling(Ts*I_day)  #number of observations is observations per day*days; is meant to keep the assessment times and the intevals equal over Ts
  dt = S_day/ I_day  #  dt in seconds;note that by using D_h=12, and dividing by I-day=24, we assume each interval has its night
  Ts_s= Ts*S_day  #transformation of season length to seconds
  ### R code rewritten
  t= seq(from=0, to=Ts, length.out = (n+1))  #times at which function is evaluated; note the extra time step for t=0; and it is now in days
  n1 = which(abs(t-Ft)==min(abs(t-Ft))) # find nearest flowering date
  #which(t == Ft)
  ### R code rewritten
  t_s= t*S_day  #transformation of time from days to seconds
  
  # prepare files and arguments
  Th_e.i = Th_e*I_day
  Th_s.i = Th_s*I_day
  Tp_s.i = Tp_s*I_day
  Tp_e.i = Tp_s.i+Tp_d*I_day # calculates the time when priming ends
  Th_d.i = Th_s.i+Th_d*I_day # calculates when defense reaches its potential
  Tco_d.i = Tco_d * I_day
  
  LAI_f= rep(0,n1)
  LAI_r= rep(0,n1)
  dLAI_tf = rep(0,n1)
  dLAI_tr = rep(0,n1)
  dLAI_mu_tf = rep(0,n1)
  dLAI_mu_tr = rep(0,n1)
  LAI_mu_f = rep(0,n1)
  LAI_mu_r = rep(0,n1)
  Pn_f= rep(0,n1)
  Pn_r= rep(0,n1)
  Pn_mu_f= rep(0,n1)
  Pn_mu_r= rep(0,n1)
  
  LAI_f[1] = Lai_f0 # start LAi values 
  LAI_r[1]= Lai_r0  # start LAi values 
  # looping
  #  for (i in 1:length(Ft)){ # do for each flowering time 
  for (j in 1:(n1-1)){  #n should thus be from t=0 to Ft (j = hours)
    if (t[j]< Ft & j < (Tp_s.i)){ # when priming is not switched on yet, the mutant will behave similar as the focal plant
      dt_out= outcomes_t(D_r=D_r,LAI_f = LAI_f[j], LAI_r = LAI_r[j], Io=Io, k=k, Lue=Lue, C=C, C_newleaf=C_newleaf,alpha=alpha, 
                         C_def_f=C_def_f,C_def_r=C_def_r,C_def_max = C_def_max,C_prim_f=C_prim_f,C_prim_r=C_prim_r,SLA=SLA,Co_r=Co_r,fPrim_f=fPrim_f,fPrim_r=fPrim_r,dt=dt, 
                         Ts=Ts_s, t=j,Th_e.i =Th_e.i, Th_s.i=Th_s.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,Def.f=Def.f,Def.r=Def.r,H.f=H.f,
                         H.r=H.r,Prim.f=Prim.f,Prim.r=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)
      Pn_f[j]  = dt_out["Pn_f"]
      Pn_r[j]  = dt_out["Pn_r"]
      LAI_f[j+1] = dt_out["LAI_f"] 
      LAI_r[j+1]= dt_out["LAI_r"]  #total= totalt-1+ time interval* function (=growth rate)
      dLAI_tf[j] = dt_out["dLAI_tf"]
      dLAI_tr[j] = dt_out["dLAI_tr"]
      LAI_mu_f[j+1] = LAI_f[j+1]
      LAI_mu_r[j+1] = LAI_r[j+1]
      Pn_mu_f[j] = Pn_f[j]
      Pn_mu_r[j] = Pn_r[j]
    } else {#(t[j]< Ft & j >= (Tp_s.i)) { # when priming is switched on the mutant will behave differently than the resident
    #  browser()
      dt_out= outcomes_t(D_r=D_r,LAI_f = LAI_f[j], LAI_r = LAI_r[j], Io=Io, k=k, Lue=Lue, C=C, C_newleaf=C_newleaf,
                         C_def_f=C_def_f,C_def_r=C_def_r,C_def_max = C_def_max,C_prim_f=C_prim_f,C_prim_r=C_prim_r,SLA=SLA,Co_r=Co_r,fPrim_f=fPrim_f,fPrim_r=fPrim_r,dt=dt, Ts=Ts_s,alpha=alpha, 
                         t=j,Th_e.i =Th_e.i, Th_s.i=Th_s.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,Def.f=Def.f,Def.r=Def.r,H.f=H.f,H.r=H.r,
                         Prim.f=Prim.f,Prim.r=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)
      Pn_f[j]  = dt_out["Pn_f"]
      Pn_r[j]  = dt_out["Pn_r"]
      LAI_f[j+1] = dt_out["LAI_f"] 
      LAI_r[j+1]= dt_out["LAI_r"]  #total= totalt-1+ time interval* function (=growth rate)
      dLAI_tf[j] = dt_out["dLAI_tf"]
      dLAI_tr[j] = dt_out["dLAI_tr"]
      # LAI production of mutant
      dt_out_mut= outcomes_t(D_r=D_r,LAI_f = LAI_mu_f[j], LAI_r = LAI_mu_r[j], Io=Io, k=k, Lue=Lue, C=C, C_newleaf=C_newleaf,alpha=alpha,
                             C_def_r = C_def_r,C_def=C_def_f_mu, C_def_max = C_def_max,C_prim_f=C_prim_f_mu,C_prim_r=C_prim_r,SLA=SLA,Co_r=Co_r_mu,fPrim_f=fPrim_f_mu,fPrim_r=fPrim_r,dt=dt, Ts=Ts_s, 
                             t=j,Th_e.i =Th_e.i, Th_s.i=Th_s.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,Def.f=Def.f_mu,Def.r=Def.r,H.f=H.f,H.r=H.r,
                             Prim.f=Prim.f_mu,Prim.r=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)
      
      LAI_mu_f[j+1] = dt_out_mut["LAI_f"]
      LAI_mu_r[j+1] = dt_out_mut["LAI_r"] 
      Pn_mu_f[j] = dt_out_mut["Pn_f"]
      Pn_mu_r[j] = dt_out_mut["Pn_r"]
      dLAI_mu_tf[j] = dt_out_mut["dLAI_tf"]
      dLAI_mu_tr[j] = dt_out_mut["dLAI_tr"]
      #print(c(Pn_mu_f[j],Pn_mu_r[j]))
    } 
  }
  #else { # seed production of focal and mutant
  seeds = Seeds_season(no=i,D_r=D_r,Pn_f=Pn_f[j],Pn_r=Pn_r[j],Th_e.i =Th_e.i, Th_s.i=Th_s.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,Ts=Ts_s,dt=dt, 
                       t=t_s,j=j,n=n,LAI_f0=LAI_f[j],LAI_r0=LAI_r[j],dLai_tf=dLAI_tf[(j)],dLai_tr = dLAI_tr[j],alpha=alpha, # t in seconds, Th_s.i in hours
                       S_day=S_day,Io=Io, k=k, Lue=Lue, C=C, C_newleaf=C_newleaf, C_def_f=C_def_f,C_def_r=C_def_r,C_def_max=C_def_max,
                       C_prim_f=C_prim_f,C_prim_r=C_prim_r,SLA=SLA,Co_r=Co_r,fPrim_f=fPrim_f,fPrim_r=fPrim_r,Def.f=Def.f,Def.r=Def.r,
                       H.f=H.f,H.r=H.r,Prim.f=Prim.f,Prim.r=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)
  seeds_mutant = Seeds_season(no=(i+1000),D_r=D_r,Pn_f=Pn_mu_f[j],Pn_r=Pn_mu_r[j],Th_e.i =Th_e.i, Th_s.i=Th_s.i,Th_d.i=Th_d.i,Tp_s.i=Tp_s.i,Tp_e.i=Tp_e.i,Tco_d.i=Tco_d.i,
                              Ts=Ts_s,dt=dt, t=t_s,j=j,n=n,LAI_f0=LAI_mu_f[j],LAI_r0=LAI_mu_r[j],dLai_tf=dLAI_tf[(j)],alpha=alpha,
                              dLai_tr = dLAI_tr[j],S_day=S_day,Io=Io, k=k, Lue=Lue, C=C, C_newleaf=C_newleaf, 
                              C_def_f=C_def_f_mu,C_def_r=C_def_r,C_def_max = C_def_max,C_prim_f=C_prim_f_mu,C_prim_r=C_prim_r,SLA=SLA,Co_r=Co_r_mu,fPrim_f=fPrim_f_mu,fPrim_r=fPrim_r,
                              Def.f=Def.f_mu,Def.r=Def.r,H.f=H.f,H.r=H.r,Prim.f=Prim.f_mu,Prim.r=Prim.r,shape=shape,responsepat=responsepat,feeding=feeding)
  
  
  # store results
  LAI_f_max = max(LAI_f)
  LAI_mu_f_max = max(LAI_mu_f)
  LAI_mu_r_max = max(LAI_mu_r)
  LAI_r_max = max(LAI_r)
  seeds_tf = seeds["seeds_tf"]
  seeds_tr = seeds["seeds_tr"]
  LAI_r_fl = seeds["LAI_r"]
  LAI_f_fl = seeds["LAI_f"]
  LAI_mu_f_fl = seeds_mutant["LAI_f"]
  seeds_mu = seeds_mutant["seeds_tf"]
  seeds_tfm = (seeds_mu-seeds_tf)/(get(mod)-get(mod_mu))
  int = (get(mod)-get(mod_mu))
  Output = data.frame(LAI_f,LAI_r,LAI_mu_f,LAI_mu_r,dLAI_tf,dLAI_tr,dLAI_mu_tf,dLAI_mu_tr,Pn_f,Pn_r,Pn_mu_f,Pn_mu_r)
  write.csv(Output,paste("Output_defense_",i,"_Defense_linear_8_delay",format(Sys.time(), "%Y%m%d"),".csv",sep=""))
  t = t[1:n1]    
  return(list(data.frame("LAI_f" =LAI_f[(n1)],"LAI_mu_f"=LAI_mu_f[(n1)], "LAI_r" =LAI_r[(n1)], LAI_f_fl,LAI_mu_f_fl, LAI_r_fl,LAI_mu_f_max,LAI_mu_r_max, LAI_f_max,LAI_r_max, seeds_tr, seeds_tf,seeds_mu,seeds_tfm,int,mod_mu,row.names = NULL ),data.frame(LAI_f,LAI_r,Pn_f,Pn_r,t))) #return both the time at which function is evaluated plus the corresponding outcomes
}

find_slope0 = function(t, slope){
  ### R code rewritten
  slope  = slope[order(t)]
  t = t[order(t)]
  
  x= slope[c(-length(slope))]	# slice without last
  z= slope[c(-1)]	#slice without first
  x3= x*z             # pos*pos and neg*neg become positive,only transition from neg to pos and reverse will lead to a negative value
  a_temp= which(x3<=0.0) #finds negative values; indicate transition
  ab=numeric()
  # R code modified 
  for (i in a_temp){  #removes divergence points (where negative becomes positive, and both an increase and decrease in variable leads to increased P
    if (slope[i]> 0.0){  #therefore, only take points where P is positive first
      # original code slope[i]>= 0.0
      ab = c(ab,t[i],">")} else {
        ab = c(ab,t[i],"<")
      }
  }
  if (length(ab)==0){ab = c(NA,NA)}
  return(ab)
}


## Function runs the model, find optimum of seed production of the community, and find game theoretical solution op ESS flowering time.   
## initiate run by setting t= 0 parameters
find_endpoints_defense = function(comb,Ft){
  print("did you check the cost function of LAI_delta and Pnet??")
  print("did you check feeding rate function?")
  for (i in 1:nrow(comb)){ # do checks before model proceeds.
    if (comb$Tp_s[i] > comb$Th_s[i]){stop("Tp_s > Th_s")} 
    if (comb$Th_e[i] < comb$Th_s[i]){stop("Th_e < Th_s")} 
    if (comb$Th_e[i] == comb$Th_s[i]){stop("Th_e is equal to Th_s")}
    if (comb$Th_e[i] > comb$Ts[i]  | comb$Th_s[i] > comb$Ts[i]){stop("Th_e or Th_s larger than season length")}
    if (comb$Ft[i] > comb$Ts[i]) stop("Flowering time larger than season length: Ft > Ts")  
    if (is.null(comb$Tco_d[i])) stop("Tco_d not defined")
    if (comb$fPrim_f[i] > 0 & !(comb$responsepat[i]=="St")){stop("fPrim > 0 and responsepattern not Faster")}
  }
  
  if (is.null(Ft)){print("Ft not specified")}
  # insert progressbar
  if (nrow(comb) > 1){
    pb <- txtProgressBar(min=1,max=nrow(comb))
  }
  
  
  out = as.data.frame(matrix(NA,nrow=nrow(comb),ncol=13))
  names(out) = c("ref","seeds_tr","seeds_tf","seeds_mu","seeds_tfm","LAI_f","LAI_r","LAI_mu_f","LAI_fl","LAI_mu_fl","LAI_f_max","LAI_r_max","LAI_mu_f_max")  
  
  for (i in 1:nrow(comb)){ # calculates the seed production and LAI over time with a given combination of parameters and stores results in out 
    
    D_f= 1  #note that we alsways take the focal plant density as 1; if one wants to model uneven ratiosn, the script has to be ammended accordingly
    D_r= comb$D[i]-D_f  ## hence, rest of the vegation is density- 1
    LAI_f0= Start_LAI(comb$S_m[i], D_f, comb$F_la[i], comb$SLA[i])  #set start LAi of focal plant
    LAI_r0= Start_LAI(comb$S_m[i], D_r, comb$F_la[i], comb$SLA[i])  #set start LAI of rest of the vegation
    # calculates the seed production and LAI over time with a given combination of parameters
    Euler_outcome = Explicit_Euler_defense(i=comb$no[i],D_r=D_r,Lai_f0=LAI_f0, Lai_r0=LAI_r0, Ts=comb$Ts[i],S_m=comb$S_m[i],F_la=comb$F_la[i],Ft=comb$Ft[i], I_day=comb$I_day[i], S_day=comb$S_day[i], 
                                           Io=comb$Io[i], k=comb$k[i], Lue=comb$Lue[i],SLA=comb$SLA[i], C=comb$C[i], C_newleaf=comb$C_newleaf[i],alpha= comb$alpha[i],
                                           C_def_r = comb$C_def_r[i],C_def_f=comb$C_def_f[i],C_def_max=comb$C_def_max[i],C_prim_f=comb$C_prim_f[i],C_prim_r=comb$C_prim_r[i],Co_r=comb$Co_r[i],fPrim_f = comb$fPrim_f[i],fPrim_r=comb$fPrim_r[i],Th_s=comb$Th_s[i],Th_e=comb$Th_e[i],Tco_d = comb$Tco_d[i],
                                           Tp_s=comb$Tp_s[i],Th_d= comb$Th_d[i],Tp_d = comb$Tp_d[i],Def.f=comb$Def.f[i],Def.r=comb$Def.r[i],H.f=comb$H.f[i],H.r=comb$H.r[i],
                                           Prim.f=comb$Prim.f[i],Prim.r=comb$Prim.r[i],mod=comb$mod[i],mut.val=comb$mut.val[i],shape=comb$shape[i],responsepat = comb$responsepat[i],feeding=comb$feeding[i])
    out[i,1] = comb$ref[i] 
    out[i,2] = max(Euler_outcome[[1]]$seeds_tr)  # seed production of resident
    out[i,3] = max(Euler_outcome[[1]]$seeds_tf)  # seed production of focal
    out[i,4] = max(Euler_outcome[[1]]$seeds_mu)  # seed production of mutant
    out[i,5] = max(Euler_outcome[[1]]$seeds_tfm) # seed df mutant and focal
    out[i,6] = Euler_outcome[[1]]$LAI_f # LAI of focal plant at end of growing season
    out[i,7] = Euler_outcome[[1]]$LAI_r # LAI of focal plant at end of growing season
    out[i,8] = Euler_outcome[[1]]$LAI_mu_f # LAI of mutant plant at end of growing season
    out[i,9] = Euler_outcome[[1]]$LAI_f_fl # LAI at the end of flowering stage of focal
    out[i,10] = Euler_outcome[[1]]$LAI_mu_f_fl # LAI at the end of flowering stage of mutant
    out[i,11] = Euler_outcome[[1]]$LAI_f_max # maximum LAI of focal plant 
    out[i,12] = Euler_outcome[[1]]$LAI_r_max # maximum LAI of resident population
    out[i,13] = Euler_outcome[[1]]$LAI_mu_f_max # maximum LAI of mutant
    # update progressbar
    #setTxtProgressBar(pb, i)  
    print(i)
  }
  # write intermediate results
  #write.csv(data.frame(comb,out),paste(format(Sys.time(), "%Y%m%d"),"_output_",eval(comb$mod[i]),".csv",sep=""))
  
  return(out)
}

find_comb = function(comb,Ft){
  # find Opt and Gt for a the combinations of parameters given the parameter of interest (comb$mod)
  # not when comb$mod = Prim.f because only true false
  out = comb[,c(40:ncol(comb))]
  comb = comb[,c(1:39)]
  comb$mod = as.character(comb$mod)
  
  if (eval(comb$mod[1]) != "Prim.f"){               
    #which.pars = which(lapply((apply(unique(comb[,-which(names(comb) == eval(comb$mod[1]))]),2,unique)),length) > 1) # retrieve the parameters that are varied along with mod
    which.pars = which(lapply((apply(unique(comb[,-grep(eval(comb$mod[1]),names(comb))]),2,unique)),length) > 1) # retrieve the parameters that are varied along with mod
    print(which.pars)
    if (Ft){which.pars = which.pars[-which(names(which.pars) =="Ft" )]} 
    print(which.pars)
    pars.comb = unique(comb[,names(which.pars)]) # extract the different combinations of parameters
    # if multiple parameters are varied
    if (length(which.pars) > 1){
      
      out.Opt = as.data.frame(matrix(NA,nrow=nrow(pars.comb),ncol=ncol(pars.comb)+3))  
      
      names(out.Opt)[1:length(which.pars)] = names(which.pars)
      names(out.Opt)[(length(which.pars)+1):(length(which.pars)+3)] = c("Opt","Gt","slopeGt")
      for (i in 1:nrow(pars.comb)){
        # take subset from out with combination of parameter values
        vec = !logical(length= nrow(comb))
        for (j in 1:ncol(pars.comb)){
          vec = vec* comb[,names(which.pars[j])] == pars.comb[i,j]
        }
        
        out.sub = out[vec,]
        par.sub = comb[vec,]
        
        max.index = which.max(out.sub$seeds_tf)
        Opt = par.sub[max.index,paste(eval(par.sub$mod[j]),sep="")] #Opt
        Gt = find_slope0(fun.name(par.sub,par.sub$mod[j]), out.sub$seeds_tfm)
        out.Opt[i,1:ncol(pars.comb)] = pars.comb[i,] # Opt
        out.Opt[i,(ncol(pars.comb)+1)] = Opt # Opt
        out.Opt[i,(ncol(pars.comb)+2):(ncol(pars.comb)+3)] = Gt  # Gt
        
      }
      # if a single parameter is varied
    } else {
      out.Opt = as.data.frame(matrix(NA,nrow=length(pars.comb),ncol=1+3))  
      names(out.Opt)[1:length(which.pars)] = names(which.pars) # maybe wrong for one parameter
      names(out.Opt)[(length(which.pars)+1):(length(which.pars)+3)] = c("Opt","Gt","slopeGt")
      for (i in 1:length(pars.comb)){
        # take subset from out with combination of parameter values
        vec = !logical(length= nrow(comb))
        vec = vec* comb[,names(which.pars[1])] == pars.comb[i]
        out.sub = out[vec,]
        par.sub = comb[vec,]
        
        max.index = which.max(out.sub$seeds_tf)
        Opt = par.sub[max.index,paste(eval(par.sub$mod[i]),sep="")] #Opt
        Gt = find_slope0(fun.name(par.sub,par.sub$mod[i]), out.sub$seeds_tfm)
        out.Opt[i,1:ncol(pars.comb)] = pars.comb[i,] # Opt
        out.Opt[i,(ncol(pars.comb)+1)] = Opt # Opt
        out.Opt[i,(ncol(pars.comb)+2):(ncol(pars.comb)+3)] = Gt  # Gt
      }
    }
  } else (out.Opt = NA) # when comb$mod[i] = Prim.f
  #    print(max(Euler_outcome$LAI_f))
  return(out.Opt)
}

# comb <- function(x, ...) {
#   lapply(seq_along(x),
#          function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
# }
# 
# # define cluster
# # cl <- makeCluster(12)
# # registerDoParallel(cl)
# 
# # do computation
# # paralllel.sim = function(x){
# # oper <- foreach(e=iter(pars, by='row'), .combine='comb', .multicombine=TRUE,
# #                 .init=list(list(), list())) %dopar% {
# #                   combine.parallel(e)
# #                 }
# #   return(oper)
# # }
# 
# # merge list 
# Merge_R_List <- function (list,ncol) {
#   a = data.frame(matrix(NA, nrow=length(list),ncol=ncol))
#   for(i in 1:length(list))    {
#     a[i,1:length(list[[i]])] = list[[i]]
#   }
#   return(a)
# }
# 
# out1 = Merge_R_List(oper[[1]],3)
# 
# 
# write.csv(out1,"d11_out_OptGt.csv")
# 
# 
# a = list.files(path = ".", pattern = "input_pars")
# read.csv(a[1])
# # 
#  LAI_delta  = cmpfun(LAI_delta)
#  Explicit_Euler_defense = cmpfun(Explicit_Euler_defense)
#  Seeds_season= cmpfun(Seeds_season)
#  Seeds_m = cmpfun(Seeds_m)
#  Seeds_season_simple =cmpfun(Seeds_season_simple)
#  find_endpoints_defense = cmpfun(find_endpoints_defense)
#  find_slope0 = cmpfun(find_slope0)
#  outcomes_t = cmpfun(outcomes_t)
#  Pnet = cmpfun(Pnet)
#  Pgross =cmpfun(Pgross)
# # 
#  
#  
