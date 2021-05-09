## ==================================================================== ##
##                                                                      ##
##                The nine-Box model program - version 1.0.0            ##
##                                                                      ##
##   ****************  Updated by YQ- on 23/10/2020  ****************   ##
##                                                                      ##
## ==================================================================== ##

## ==================================================================== ##
setwd("E:/UoB PhD/Second year_phD/20200620_9box_model/code/16box_v1.0.0_RCS")

## ************************* library packages ************************* ##
library(tidyverse)
## library(deSolve)

## ************************* remove history *************************** ##
rm(list = ls(all = T))

## =======================================================================
##  *********************** MODEL CONFIGURATION *********************** ##
## =======================================================================

## ******************* STREET (BOX) configuration ********************* ##
## input:
## 	 nbv,nbh num of boxes on verical and horizontal dir. 
## 	 nbox num of boxes in total
##	 h1,h2 the height of the boxes
##	 l1,l2 the width  of the boxes
nbv = 16; nbh = 4; nbox = 16;

Htotal = 36; Ltotal = 18;

alphaH = array(c(0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625),c(nbv,1));
betaL = array(c(0.25,0.25,0.25,0.25),c(nbh,1)) # The summation is equal to 1.0!

#
h = array(0,c(nbv,1)); l = array(0,c(nbh,1)) # The summation is equal to 1.0!
h[1] = alphaH[1]*Htotal; h[2] = alphaH[2]*Htotal; h[3] = alphaH[3]*Htotal; h[4] = alphaH[4]*Htotal
l[1] = betaL[1]*Ltotal; l[2] = betaL[2]*Ltotal; l[3] = betaL[3]*Ltotal; l[4] = betaL[4]*Ltotal

## ************************* Loop parameters ************************** ##
## input:
## 	 Numspe total num of species + passive scalar (here is 52th); NoS num of handling species
## 	 deltaTime,idtlong slow chemistry timestep; dtchem fast chemistry timestep
##	 counter,nshort control fast reactions
##	 nprint outputs every min = nprint*timestep
##   Carr output matrix of Numspe.th species in box[nbv,nbh]
##   nosteps total timestep (total time = nosteps * deltaTime sec > 1800s)

Numspe = 52; NoS = 1; Carr = array(0,c(nbv,nbh,Numspe));

deltaTime = 0.03; idtlong = deltaTime; dtchem = 0.003; 

nprint = 2000; nosteps = 420000;
# test; test=(int)min(10.0,max(2.38,1.0)); iDTsno3; iNsno3; ino3counter = 0

## =======================================================================
##  ************************* Input Document ************************** ##
## =======================================================================
## input:
## 	 cB background concentrations on the rooftop
## 	 trafficin traffic amounts at leeward [1,1] and windward [1,2]
##   cloudin cloud cover
##   beta photolysis rates
##   PS=PS+EmiPS*deltaTime;
## PS=PS+EmiPS*deltaTime*traScale; # what is trascale?????
cBgarr = array(0,c(Numspe,1))
beta = 1.0;
gammaEmi = array(c(0.25,0.25,0.25,0.25),c(nbh,1));
trafficin = array(c(1500,1500,1500,1500,1500,1500,1500,1500,1500,1500, 
                    1500,1500,1500,1500,1500,1500,1500,1500,1500,1500),c(10,nbh))
trafficin = t(t(trafficin)*c(gammaEmi[1],gammaEmi[2],gammaEmi[3],gammaEmi[4]))
#
wtin = array(c(0.021,0.021,0.021,0.021,0.021,0.021,0.021,0.021,0.021,0.021),c(10,1)); 
Scalewt = 1.0; wtin=Scalewt*wtin;
#
wd = 0; cloudin = array(0,c(10));

##  ********************** Control vortex direction *********************** ##
uein = array(0,c(nbv+1,nbh))	# horizontal exchange vel
wein = array(0,c(nbv+1,nbh))	# vertical exchange vel
uain = array(0,c(nbv+1,nbh))  # horizontal adv vel
wain = array(0,c(nbv+1,nbh))	# vertical adv vel
if (wd > 90 & wd <= 270){
    wain[2,4] = 0.1468; wain[3,4] = 0.2537; wain[4,4] = 0.1805;
    wain[2,3] = 0.0938; wain[3,3] = 0.1850; wain[4,3] = 0.1335;
    wain[2,2] =-0.0712; wain[3,2] =-0.1252; wain[4,2] =-0.0642;
    wain[2,1] =-0.1670; wain[3,1] =-0.2888; wain[4,1] =-0.2597;
    wein[2,4] = 0.0111; wein[3,4] = 0.0119; wein[4,4] = 0.0174; wein[5,4] = 0.0097;
    wein[2,3] = 0.0136; wein[3,3] = 0.0041; wein[4,3] = 0.0194; wein[5,3] = 0.0224;
    wein[2,2] = 0.0110; wein[3,2] = 0.0115; wein[4,2] = 0.0066; wein[5,2] = 0.0228;
    wein[2,1] = 0.0815; wein[3,1] = 0.1691; wein[4,1] = 0.0649; wein[5,1] = 0.0358;
    uain[1,4] = 0.1684; uain[1,3] = 0.2872; uain[1,4] = 0.2171;
    uain[2,4] = 0.0952; uain[2,3] = 0.1727; uain[2,4] = 0.1430;
    uain[3,4] =-0.0753; uain[3,3] =-0.1304; uain[3,4] =-0.0672;
    uain[4,4] =-0.1751; uain[4,3] =-0.2812; uain[4,4] =-0.2111;
    uein[1,4] = 0.0216; uein[1,3] = 0.0102; uein[1,4] = 0.0167;
    uein[2,4] = 0.0180; uein[2,3] = 0.0233; uein[2,4] = 0.0107;
    uein[3,4] = 0.0200; uein[3,3] = 0.0588; uein[3,4] = 0.5195;
    uein[4,4] = 0.0174; uein[4,3] = 0.0400; uein[4,4] = 0.0021;
  # counterclockwise
  } else {
    wain[2,1] = 0.1468; wain[3,1] = 0.2537; wain[4,1] = 0.1805;
    wain[2,2] = 0.0938; wain[3,2] = 0.1850; wain[4,2] = 0.1335;
    wain[2,3] =-0.0712; wain[3,3] =-0.1252; wain[4,3] =-0.0642;
    wain[2,4] =-0.1670; wain[3,4] =-0.2888; wain[4,4] =-0.2597;
    wein[2,1] = 0.0111; wein[3,1] = 0.0119; wein[4,1] = 0.0174; wein[5,1] = 0.0097;
    wein[2,2] = 0.0136; wein[3,2] = 0.0041; wein[4,2] = 0.0194; wein[5,2] = 0.0224;
    wein[2,3] = 0.0110; wein[3,3] = 0.0115; wein[4,3] = 0.0066; wein[5,3] = 0.0228;
    wein[2,4] = 0.0815; wein[3,4] = 0.1691; wein[4,4] = 0.0649; wein[5,4] = 0.0358;
    uain[1,2] =-0.1684; uain[1,3] =-0.2872; uain[1,4] =-0.2171;
    uain[2,2] =-0.0952; uain[2,3] =-0.1727; uain[2,4] =-0.1430;
    uain[3,2] = 0.0753; uain[3,3] = 0.1304; uain[3,4] = 0.0672;
    uain[4,2] = 0.1751; uain[4,3] = 0.2812; uain[4,4] = 0.2111;
    uein[1,2] = 0.0216; uein[1,3] = 0.0102; uein[1,4] = 0.0167;
    uein[2,2] = 0.0180; uein[2,3] = 0.0233; uein[2,4] = 0.0107;
    uein[3,2] = 0.0200; uein[3,3] = 0.0588; uein[3,4] = 0.5195;
    uein[4,2] = 0.0174; uein[4,3] = 0.0400; uein[4,4] = 0.0021;
  # clockwise
  }

##  *********************************************************************** ##
# note: u=ue; w=we; U=ua; W=wa
# T0 = 263; # Temperature in Kelvin
# Q0 = 1000; # Solar radiation in W/m2

## =======================================================================
##  ********************* Emission factor control ********************* ##
## =======================================================================
## input:
## 	 Emiarr[conc,1/2] vehicle emission factors at 01 or 02 rd, 52th refers to passive scalar
Emiarr = array(0,c(Numspe,nbh))
#
EmiNOx = 1000;
EmiVOC = 791;
EmiCO = 3593*0.3*0.3/(h[1]*l[1]); # equal to 337 ug/ms
EmiPS = 1000*0.3*0.3/(h[1]*l[1]);

alphaNO2 = 0.1;
ScaleEmiNOx = 0.1; # ScaleEmiNOx
noScaleEmiNOx = 10; # noScaleEmiNOx
ScaleEmiVOC = 0.1; # ScaleEmiVOC
noScaleEmiVOC = 10; # noScaleEmiVOC

EmiNOx=EmiNOx*0.3*0.3/(h[1]*l[1]);
EmiVOC=EmiVOC*0.3*0.3/(h[1]*l[1]);
EmiNOx=noScaleEmiNOx*ScaleEmiNOx*EmiNOx;
EmiVOC=noScaleEmiVOC*ScaleEmiVOC*EmiVOC;

Emiarr[1,] = EmiNOx*(1.0-alphaNO2); # EmiNO
Emiarr[2,] = EmiNOx*alphaNO2; # EmiNO2
Emiarr[8,] = EmiNOx*0.008; # EmiHONO   REF:kirchstetter, Kurtenbach...
Emiarr[6,] = EmiCO; # EmiCO
Emiarr[12,] = (347.0/791.0)*EmiVOC; # EmiC2H4
Emiarr[13,] = (150.0/791.0)*EmiVOC; # EmiC3H6
Emiarr[23,] = (98.0/791.0)*EmiVOC; # EmiCH3CHO
Emiarr[26,] = (196.0/791.0)*EmiVOC; # EmiHCHO
Emiarr[52,] = EmiPS; # EmiHCHO

## =======================================================================
##  ***********************  Initial condition  *********************** ## 
## =======================================================================

# Initial species NO.1-10
iNO = array(1.33,c(nbv,nbh));
iNO2 = array(6.78,c(nbv,nbh));
iNO3 = array(0,c(nbv,nbh));
iO3 = array(49.99,c(nbv,nbh));
iH2 = array(0.00,c(nbv,nbh));

iCO = array(119.70,c(nbv,nbh));
iH2O2 = array(0.00,c(nbv,nbh));
iHONO  = array(0.00,c(nbv,nbh));
iHNO3 = array(2.00,c(nbv,nbh));
iHO2NO2 = array(0.00,c(nbv,nbh));

# Initial species NO.11-20
iCH4 = array(1800.00,c(nbv,nbh));
iC2H4 = array(0.91,c(nbv,nbh));
iC3H6 = array(0.29,c(nbv,nbh));
iC5H8 = array(0.28,c(nbv,nbh));
iC2H5OH = array(2.37,c(nbv,nbh));

iCH3NO3 = array(0.00,c(nbv,nbh));
iOH = array(0.00,c(nbv,nbh));
iHO2 = array(0.00,c(nbv,nbh));
iCH3O2 = array(0.00,c(nbv,nbh));
iCH3OH = array(7.38,c(nbv,nbh));

# Initial species NO.21-30
iHOCH2CH2O2 = array(0.00,c(nbv,nbh));
iCH3CO2H = array(0.00,c(nbv,nbh));
iCH3CHO = array(2.98,c(nbv,nbh));
iRU14O2 = array(0.00,c(nbv,nbh));
iUCARB10 = array(0.00,c(nbv,nbh));

iHCHO = array(3.14,c(nbv,nbh));
iCH3CO3 = array(0.00,c(nbv,nbh));
iHOCH2CHO = array(0.00,c(nbv,nbh));
iRN9O2 = array(0.00,c(nbv,nbh));
iHCOOH = array(0.00,c(nbv,nbh));

# Initial species NO.31-40
iCARB6 = array(0.00,c(nbv,nbh));
iUCARB12 = array(0.00,c(nbv,nbh));
iRU12O2 = array(0.00,c(nbv,nbh));
iCARB7 = array(0.00,c(nbv,nbh));
iRU10O2 = array(0.00,c(nbv,nbh));

iHOC2H4NO3 = array(0.00,c(nbv,nbh));
iRN9NO3 = array(0.00,c(nbv,nbh));
iRU14NO3 = array(0.00,c(nbv,nbh));
iCH3OOH = array(0.00,c(nbv,nbh));
iHOC2H4OOH = array(0.00,c(nbv,nbh));

# Initial species NO.41-51
iRN9OOH = array(0.00,c(nbv,nbh));
iCH3CO3H = array(0.00,c(nbv,nbh));
iHOCH2CO3H = array(0.00,c(nbv,nbh));
iRU14OOH = array(0.00,c(nbv,nbh));
iRU12OOH = array(0.00,c(nbv,nbh));

iRU10OOH = array(0.00,c(nbv,nbh));
iHOCH2CO3 = array(0.00,c(nbv,nbh));
iPAN = array(0.46,c(nbv,nbh));
iPHAN = array(0.00,c(nbv,nbh));
iRU12PAN = array(0.00,c(nbv,nbh));

iMPAN = array(0.00,c(nbv,nbh))

# Initial Passive scalar NO.52
iPS = array(0.00,c(nbv,nbh))

##  ****************  initial concentrations for input **************** ## 
Carr[,,1] = iNO;
Carr[,,2] = iNO2;
Carr[,,3] = iNO3;
Carr[,,4] = iO3;
Carr[,,5] = iH2;

Carr[,,6] = iCO;
Carr[,,7] = iH2O2;
Carr[,,8] = iHONO;
Carr[,,9] = iHNO3;
Carr[,,10] = iHO2NO2;

Carr[,,11] = iCH4;
Carr[,,12] = iC2H4;
Carr[,,13] = iC3H6;
Carr[,,14] = iC5H8;
Carr[,,15] = iC2H5OH;

Carr[,,16] = iCH3NO3;
Carr[,,17] = iOH;
Carr[,,18] = iHO2;
Carr[,,19] = iCH3O2;
Carr[,,20] = iCH3OH;

Carr[,,21] = iHOCH2CH2O2;
Carr[,,22] = iCH3CO2H;
Carr[,,23] = iCH3CHO;
Carr[,,24] = iRU14O2;
Carr[,,25] = iUCARB10;

Carr[,,26] = iHCHO;
Carr[,,27] = iCH3CO3;
Carr[,,28] = iHOCH2CHO;
Carr[,,29] = iRN9O2;
Carr[,,30] = iHCOOH;

Carr[,,31] = iCARB6;
Carr[,,32] = iUCARB12;
Carr[,,33] = iRU12O2;
Carr[,,34] = iCARB7;
Carr[,,35] = iRU10O2;

Carr[,,36] = iHOC2H4NO3;
Carr[,,37] = iRN9NO3;
Carr[,,38] = iRU14NO3;
Carr[,,39] = iCH3OOH;
Carr[,,40] = iHOC2H4OOH;

Carr[,,41] = iRN9OOH;
Carr[,,42] = iCH3CO3H;
Carr[,,43] = iHOCH2CO3H;
Carr[,,44] = iRU14OOH;
Carr[,,45] = iRU12OOH;

Carr[,,46] = iRU10OOH;
Carr[,,47] = iHOCH2CO3;
Carr[,,48] = iPAN;
Carr[,,49] = iPHAN;
Carr[,,50] = iRU12PAN;

Carr[,,51] = iMPAN

Carr[,,52] = iPS
## **************************  COMPLETE!  ***************************** ##

## =======================================================================
## ******* load chemistry & exchange module & timescale module ******** ##
## =======================================================================
  source("16Exchange_module_v1.0.0.R")
  source("16Chemi_module_RCS_v1.0.0.R")
## =======================================================================
## ************************** The time loop *************************** ##
## created by YQ- on 26/10/2020                                         ##
## =======================================================================
for (nint in 1:nosteps)
{  
  ## Setup initial/background concentrations with 30 min spin-up
  if (nint <= 1800.00/deltaTime) 
    {
    ## Chemistry
      ## NO.1 beta =1.0 for spin-up
      Carr = chemi(Carr, deltaTime, beta = 1.0, dtchem)
      for (NoS in 1:Numspe)
      {
        cBgarr[NoS] = Carr[1,1,NoS]
      }
      ## Calculate chemical timescale
      ## Tchem = chemti(Carr,beta)
      ## nint = nint +1
      
    } else
    { ## Adding emissions & exchange
      # Note: for one hr run, the nint step is 3600/0.03 (deltaTime) = 120000!
      myquot = nint %/% 120000;
      myrem = nint %% 120000;
      ## update time varying ue;we;ua;wa 
      wt = wtin[myquot+1]+(myrem/120000)*(wtin[myquot+2]-wtin[myquot+1]);
      ue = uein*(wt/0.021)
      we = wein*(wt/0.021)
      ua = uain*(wt/0.021)
      wa = wain*(wt/0.021)
      ## update time varying emissions
      traffic=trafficin[myquot+1,]+(myrem/120000)*(trafficin[myquot+2,]-trafficin[myquot+1,]);
      traScale=traffic/1500.0; # Emission profile refers to 1500 vehicles per hour
      ## update time varying cloud cover
      cloud=cloudin[myquot+1]+(myrem/120000)*(cloudin[myquot+2]-cloudin[myquot+1]);
      beta=(8.0-cloud)/8.0; # update photolysis rate
 
      ## ************** update the concentration of Carr **************** ##
      ## Chemistry
      Carr = chemi(Carr, deltaTime, beta, dtchem)
      
      ## **************** add emissions into the model ****************** ##
      ## 
      for (NoS in 1:Numspe)
      {
        Carr[1,1,NoS] = Carr[1,1,NoS]+Emiarr[NoS,1]*deltaTime*traScale[1]; # 1st road!
        Carr[1,2,NoS] = Carr[1,2,NoS]+Emiarr[NoS,2]*deltaTime*traScale[2]; # 2nd road!
        Carr[1,3,NoS] = Carr[1,3,NoS]+Emiarr[NoS,3]*deltaTime*traScale[3]; # 3rd road!
        Carr[1,4,NoS] = Carr[1,4,NoS]+Emiarr[NoS,4]*deltaTime*traScale[4]; # 4th road!
      } 
      
      ## ******** add exchanges with background into the model ********** ##
      Carr = exclock(u = ue, w = we, U = ua, W = wa, Carr, cBgarr,
               h,l, nbox, nbv, nbh,
               deltaTime, NoS, Numspe)
      }

      ## *******************  COMPLETE one round!  ********************** ##
  if (nint %% nprint == 0)
  {
      write(Carr, file = "16box_RCSs_BASE_yzywt.csv", ncolumns = if(is.character(Carr)) 1 else 832, 
            append = T, sep = ",")
  }
}
## *************************** !!!COMPLETE!!! *************************** ##
