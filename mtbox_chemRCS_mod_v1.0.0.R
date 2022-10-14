## =======================================================================
##  *********** load chemistry module for multi-box model  ***********  #
## Created by YQ- on 23/07/2020 based on the work of 
## Zhong, J, Cai, X & Bloss, W 2017, Large eddy simulation of reactive pollutants in a deep urban street canyon: Coupling dynamics with O3-NOx-VOC chemistry, Environmental Pollution. 224: 171-184.

## =======================================================================
chemi <- function (Carr, deltaTime, beta, dtchem) 
{
  idtlong = deltaTime;
  counter = 0;
  nshort = idtlong/dtchem;
  ## =======================================================================
  ##  ********* declaration inital concentration for chemistry *********  ##
  ##  updated by YQ- on 23/08/2020 
  ## =======================================================================
  iNO = Carr[,,1]; 
  iNO2 = Carr[,,2]; 
  iNO3 = Carr[,,3]; 
  iO3 = Carr[,,4];
  iH2 = Carr[,,5];
  
  iCO = Carr[,,6];
  iH2O2 = Carr[,,7];
  iHONO = Carr[,,8];
  iHNO3 = Carr[,,9];
  iHO2NO2 = Carr[,,10];
  
  iCH4 = Carr[,,11];
  iC2H4 = Carr[,,12];
  iC3H6 = Carr[,,13];
  iC5H8 = Carr[,,14];
  iC2H5OH = Carr[,,15];
  
  iCH3NO3 = Carr[,,16];
  iOH =	Carr[,,17];
  iHO2 = Carr[,,18];
  iCH3O2 = Carr[,,19];
  iCH3OH = Carr[,,20];
  
  iHOCH2CH2O2 =	Carr[,,21];
  iCH3CO2H = Carr[,,22];
  iCH3CHO =	Carr[,,23];
  iRU14O2 = Carr[,,24];
  iUCARB10 = Carr[,,25];
  
  iHCHO = Carr[,,26];
  iCH3CO3 = Carr[,,27];
  iHOCH2CHO = Carr[,,28];
  iRN9O2 = Carr[,,29];
  iHCOOH = Carr[,,30];
  
  iCARB6 = Carr[,,31];
  iUCARB12 = Carr[,,32];
  iRU12O2 = Carr[,,33];
  iCARB7 = Carr[,,34];
  iRU10O2 = Carr[,,35];
  
  iHOC2H4NO3 = Carr[,,36];
  iRN9NO3 = Carr[,,37];
  iRU14NO3 = Carr[,,38];
  iCH3OOH = Carr[,,39];
  iHOC2H4OOH = Carr[,,40];
  
  iRN9OOH =	Carr[,,41];
  iCH3CO3H = Carr[,,42];
  iHOCH2CO3H = Carr[,,43];
  iRU14OOH = Carr[,,44];
  iRU12OOH = Carr[,,45];
  
  iRU10OOH = Carr[,,46];
  iHOCH2CO3 = Carr[,,47];
  iPAN = Carr[,,48];
  iPHAN = Carr[,,49];
  iRU12PAN = Carr[,,50];
  
  iMPAN = Carr[,,51];
  
  ## =======================================================================
  ##  *************** declaration for reaction constants **************** ##
  ##  updated by YQ- on 23/07/2020
  ## =======================================================================
  ick = array(rep(c(3.40e-6,	# ick1
                4.01e-4,	# ick2
                2.63e-9,	# ick3
                6.56e-1,	# ick4
                1.72e-3,	# ick5
                1.49e-4,	# ick6
                5.06e-3,	# ick7
                4.21e-2,	# ick8
                4.86e-5,	# ick9
                2.82e+0,	# ick10
                8.74e-2,	# ick11
                6.92e-2,	# ick12
                2.54e-1,	# ick13
                3.08e-1,	# ick14
                5.01e-1,	# ick15
                2.27e-1,	# ick16
                3.59e-2,	# ick17
                3.74e-2,	# ick18
                1.20e-1,	# ick19
                2.58e-2,	# ick20
                4.08e-3,	# ick21
                7.11e-6,	# ick22
                9.20e-3,	# ick23
                2.34e-2,	# ick24
                1.83e-1,	# ick25
                2.02e-3,	# ick26
                6.30e-7,	# ick27
                1.39e-4,	# ick28
                2.00e-1,	# ick29
                7.19e-1,	# ick30
                4.46e-9,	# ick31
                2.99e-8, 	# ick32
                8.18e-8,	# ick33
                1.45e-7,	# ick34
                2.58e+0,	# ick35
                7.76e-8,	# ick36
                2.10e-7,	# ick37
                3.05e-5,	# ick38
                4.61e-5,	# ick39
                5.07e-6, 	# ick40
                2.35e-1,	# ick41
                4.02e-1,	# ick42
                2.31e-2,	# ick43
                7.24e-2,	# ick44
                9.23e-3,	# ick45
                1.13e-2,	# ick46
                2.00e-2,	# ick47
                1.95e-1, 	# ick48
                1.68e-1,	# ick49
                4.84e-2,	# ick50
                2.13e-1,	# ick51
                5.10e-1,	# ick52
                5.10e-1,	# ick53
                4.93e-2,	# ick54
                1.46e-1,	# ick55
                1.52e-1,	# ick56
                6.52e-2,	# ick57
                1.09e-1,	# ick58
                6.52e-2,	# ick59
                4.35e-2,	# ick60
                1.95e-4,	# ick61
                1.09e-3,	# ick62
                4.56e-3,	# ick63
                2.17e-2, 	# ick64
                1.52e-1,	# ick65
                3.62e-1,	# ick66
                3.20e-1,	# ick67
                3.75e-1,	# ick68
                3.75e-1,	# ick69
                4.74e-1,	# ick70
                4.35e-1,	# ick71
                3.85e-1, 	# ick72
                6.22e-3,	# ick73
                6.32e-3,	# ick74
                6.32e-3,	# ick75
                1.12e-2,	# ick76
                2.20e-2,	# ick77
                2.50e-1,	# ick78
                2.50e-1,	# ick79
                1.08e-2,	# ick80
                3.20e-2,	# ick81
                3.51e-2,	# ick82
                1.50e-2,	# ick83
                2.50e-2,	# ick84
                1.50e-2,	# ick85
                1.00e-2,	# ick86
                3.36e-6,	# ick87
                1.77e-5,	# ick88
                1.62e-5,	# ick89
                1.26e-4,	# ick90
                1.62e-5,	# ick91
                7.51e-2,	# ick92
                6.26e-1,	# ick93
                4.21e-8,	# ick94
                2.93e-8,	# ick95
                2.50e-1,	# ick96
                4.31e-1,	# ick97
                1.13e+0,	# ick98
                5.35e-7,	# ick99
                6.61e-8,	# ick100
                8.96e-7,	# ick101
                9.33e-3,	# ick102
                2.73e-2,	# ick103
                3.28e-2,	# ick104
                1.39e+0,	# ick105
                5.44e-6,	# ick106
                5.44e-6,	# ick107
                5.44e-6,	# ick108
                1.37e-6,	# ick109
                4.07e-6,	# ick110
                5.44e-6,	# ick111
                5.44e-6,	# ick112
                5.44e-6,	# ick113
                5.44e-6,	# ick114
                9.10e-1,	# ick115
                4.79e-1, 	# ick116
                9.27e-2,	# ick117
                1.55e-1,	# ick118
                1.88e+0,	# ick119
                7.51e-1,	# ick120
                7.51e-1,	# ick121
                5.34e-1,	# ick122
                6.26e-1,	# ick123
                2.68e-1,	# ick124
                1.51e-4,	# ick125
                2.68e-1,	# ick126
                1.51e-4,	# ick127
                2.59e-3,	# ick128
                2.81e-2,	# ick129
                1.63e-2,	# ick130
                1.51e-4,	# ick131
                1.10e-2,	# ick132
                1.51e-4,	# ick133
                9.02e-2,	# ick134
                6.31e-1,	# ick135
                7.65E-07	# ick136
  ),each = nbv*nbh),c(nbv,nbh,136));

  irk = array(rep(c(6.22e-3,	# irk73
                    6.32e-3,  # irk74
                    6.32e-3,	# irk75
                    1.12e-2,	# irk76
                    2.20e-2,	# irk77
                    2.50e-1,	# irk78
                    2.50e-1,  # irk79
                    1.08e-2,	# irk80
                    3.20e-2,	# irk81
                    3.51e-2,	# irk82
                    1.50e-2,	# irk83
                    2.50e-2,	# irk84
                    1.50e-2,	# irk85
                    1.00e-2  # irk86
  ),each = nbv*nbh),c(nbv,nbh,14));
  
  rep(c(1,2,3,4,5),each=10)
  
  iRO2 = iCH3O2+iHOCH2CH2O2+iRN9O2+iCH3CO3+iHOCH2CO3+iRU14O2+iRU12O2+iRU10O2
  
  ##   ********************* RO2 permutation reactions ************************ # 
  ick[,,73]  = irk[,,1]*iRO2; ick[,,74]  = irk[,,2]*iRO2; ick[,,75]  = irk[,,3]*iRO2; ick[,,76]  = irk[,,4]*iRO2; 
  ick[,,77]  = irk[,,5]*iRO2; ick[,,78]  = irk[,,6]*iRO2; ick[,,79]  = irk[,,7]*iRO2; ick[,,80]  = irk[,,8]*iRO2;
  ick[,,81]  = irk[,,9]*iRO2; ick[,,82]  = irk[,,10]*iRO2; ick[,,83]  = irk[,,11]*iRO2; ick[,,84]  = irk[,,12]*iRO2;
  ick[,,85]  = irk[,,13]*iRO2; ick[,,86]  = irk[,,14]*iRO2;
  
  ##   ********************* calculation of photolysis rates ****************** # 
  # 0 photolysis reaction 0-0
  ick[,,1] = beta*ick[,,1];
  # 1 photolysis reaction 1-6
  ick[,,22] = beta*ick[,,22];ick[,,23]=beta*ick[,,23];ick[,,24]=beta*ick[,,24];
  ick[,,25] = beta*ick[,,25];ick[,,26]=beta*ick[,,26];ick[,,27]=beta*ick[,,27];
  # 2 photolysis reaction 7-9
  ick[,,38] = beta*ick[,,38];ick[,,39]=beta*ick[,,39];ick[,,40]=beta*ick[,,40];
  # 3 photolysis reaction 10-14
  ick[,,87] = beta*ick[,,87];ick[,,88]=beta*ick[,,88];ick[,,89]=beta*ick[,,89];
  ick[,,90] = beta*ick[,,90];ick[,,91]=beta*ick[,,91];
  # 4 photolysis reaction 15-15
  ick[,,101] = beta*ick[,,101];
  # 5 photolysis reaction 16-24
  ick[,,106]=beta*ick[,,106];ick[,,107]=beta*ick[,,107];ick[,,108]=beta*ick[,,108];
  ick[,,109]=beta*ick[,,109];ick[,,110]=beta*ick[,,110];ick[,,111]=beta*ick[,,111];
  ick[,,112]=beta*ick[,,112];ick[,,113]=beta*ick[,,113];ick[,,114]=beta*ick[,,114];
  
  ## =======================================================================
  ##  *************** calculation of loss & prod terms ****************** ##
  ##  updated by YQ- on 23/07/2020
  ## =======================================================================
  
  ## **********************Calculation of loss terms ******************** ##
  # Solver Species no.1: NO
  iLoNO = ick[,,2]*iO3+2*ick[,,3]*iNO+ick[,,4]*iNO3+ick[,,13]*iOH+ick[,,16]*iHO2+ick[,,48]*iCH3O2+ick[,,49]*iHOCH2CH2O2+ick[,,50]*iHOCH2CH2O2+
    ick[,,51]*iRN9O2+ick[,,52]*iCH3CO3+ick[,,53]*iHOCH2CO3+ick[,,54]*iRU14O2+ick[,,55]*iRU14O2+ick[,,56]*iRU12O2+ick[,,57]*iRU12O2+
    ick[,,58]*iRU10O2+ick[,,59]*iRU10O2+ick[,,60]*iRU10O2+ick[,,61]*iCH3O2+ick[,,62]*iHOCH2CH2O2+ick[,,63]*iRN9O2+
    ick[,,64]*iRU14O2;
  # Solver Species no.2: NO2
  iLoNO2 = ick[,,14]*iOH+ick[,,17]*iHO2+ick[,,136]*iO3+
    ick[,,23]+ick[,,124]*iCH3CO3+ick[,,126]*iHOCH2CO3+ick[,,130]*iRU12O2+ick[,,132]*iRU10O2;
  # Solver Species no.4: O3
  iLoO3 = ick[,,1]+ick[,,2]*iNO+ick[,,5]*iOH+ick[,,9]*iHO2+ick[,,31]*iC2H4+ick[,,32]*iC2H4+ick[,,33]*iC3H6+ick[,,34]*iC3H6+ick[,,36]*iC5H8+
    ick[,,37]*iC5H8+ick[,,94]*iUCARB10+ick[,,95]*iUCARB10+ick[,,99]*iUCARB12+ick[,,100]*iUCARB12+ick[,,136]*iNO2;
  # Solver Species no.5: H2
  iLoH2 = ick[,,6]*iOH;
  # Solver Species no.6: CO
  iLoCO = ick[,,7]*iOH;
  # Solver Species no.7: H2O2
  iLoH2O2 = ick[,,8]*iOH+ick[,,22];
  # Solver Species no.8: HONO
  iLoHONO = ick[,,20]*iOH+ick[,,26];
  # Solver Species no.9: HNO3
  iLoHNO3 = ick[,,21]*iOH+ick[,,27];
  # Solver Species no.10: HO2NO2
  iLoHO2NO2 = ick[,,18]+ick[,,19]*iOH;
  # Solver Species no.11: CH4
  iLoCH4 = ick[,,28]*iOH;
  # Solver Species no.12: C2H4
  iLoC2H4 = ick[,,29]*iOH+ick[,,31]*iO3+ick[,,32]*iO3;
  # Solver Species no.13: C3H6
  iLoC3H6 = ick[,,30]*iOH+ick[,,33]*iO3+ick[,,34]*iO3;
  # Solver Species no.14: C5H8
  iLoC5H8 = ick[,,35]*iOH+ick[,,36]*iO3+ick[,,37]*iO3;
  # Solver Species no.15: C2H5OH
  iLoC2H5OH = ick[,,44]*iOH+ick[,,45]*iOH;
  # Solver Species no.16: CH3NO3
  iLoCH3NO3 = ick[,,101]+ick[,,102]*iOH;     
  # Solver Species no.20: CH3OH
  iLoCH3OH = ick[,,43]*iOH;
  # Solver Species no.22: CH3CO2H
  iLoCH3CO2H = ick[,,47]*iOH;
  # Solver Species no.23: CH3CHO
  iLoCH3CHO = ick[,,40]+ick[,,42]*iOH;
  # Solver Species no.25: UCARB10
  iLoUCARB10 = ick[,,89]+ick[,,93]*iOH+ick[,,94]*iO3+ick[,,95]*iO3;
  # Solver Species no.26: HCHO
  iLoHCHO = ick[,,38]+ick[,,39]+ick[,,41]*iOH;
  # Solver Species no.28: HOCH2CHO
  iLoHOCH2CHO = ick[,,88]+ick[,,96]*iOH;
  # Solver Species no.30: HCOOH
  iLoHCOOH = ick[,,46]*iOH;
  # Solver Species no.31: CARB6
  iLoCARB6 = ick[,,90]+ick[,,97]*iOH;
  # Solver Species no.32: UCARB12
  iLoUCARB12 = ick[,,91]+ick[,,98]*iOH+ick[,,99]*iO3+ick[,,100]*iO3;
  # Solver Species no.34: CARB7
  iLoCARB7 = ick[,,87]+ick[,,92]*iOH;
  # Solver Species no.35: RU10O2
  #  corrected on 20/07 iLoRU10O2 = ick[,,58]*iNO+ick[,,59]*iNO+ick[,,60]*iNO+ick[,,72]*iHO2+ick[,,84]+ick[,,85]+ick[,,86]+ick[,,132]*iNO2;
  # Solver Species no.36: HOC2H4NO3
  iLoHOC2H4NO3 = ick[,,103]*iOH;
  # Solver Species no.37: RN9NO3
  iLoRN9NO3 = ick[,,104]*iOH;
  # Solver Species no.38: RU14NO3
  iLoRU14NO3 = ick[,,105]*iOH;
  # Solver Species no.39: CH3OOH
  iLoCH3OOH = ick[,,106]+ick[,,115]*iOH+ick[,,116]*iOH;
  # Solver Species no.40: HOC2H4OOH
  iLoHOC2H4OOH = ick[,,113]+ick[,,122]*iOH;
  # Solver Species no.41: RN9OOH
  iLoRN9OOH = ick[,,114]+ick[,,123]*iOH;
  # Solver Species no.42: CH3CO3H
  iLoCH3CO3H = ick[,,107]+ick[,,117]*iOH;
  # Solver Species no.43: HOCH2CO3H
  iLoHOCH2CO3H = ick[,,108]+ick[,,118]*iOH;
  # Solver Species no.44: RU14OOH
  iLoRU14OOH = ick[,,109]+ick[,,110]+ick[,,119]*iOH;
  # Solver Species no.45: RU12OOH
  iLoRU12OOH = ick[,,111]+ick[,,120]*iOH;
  # Solver Species no.46: RU10OOH
  iLoRU10OOH = ick[,,112]+ick[,,121]*iOH;
  # Solver Species no.48: PAN
  iLoPAN = ick[,,125]+ick[,,128]*iOH;
  # Solver Species no.49: PHAN
  iLoPHAN = ick[,,127]+ick[,,129]*iOH;
  # Solver Species no.50: RU12PAN
  iLoRU12PAN = ick[,,131]+ick[,,135]*iOH;
  # Solver Species no.51: MPAN
  iLoMPAN = ick[,,133]+ick[,,134]*iOH;
  
  ## ***************** Calculation of production terms ******************** ##
  
  # solver species no.1: NO
  iPoNO = ick[,,23]*iNO2+ick[,,24]*iNO3+ick[,,26]*iHONO;
  # solver species no.2: NO2
  iPoNO2 = ick[,,2]*iNO*iO3+2*ick[,,3]*iNO*iNO+2*ick[,,4]*iNO*iNO3+ick[,,15]*iOH*iNO3+ick[,,16]*iHO2*iNO+ick[,,18]*iHO2NO2+
    ick[,,19]*iOH*iHO2NO2+ick[,,20]*iOH*iHONO+ick[,,25]*iNO3+ick[,,27]*iHNO3+ick[,,48]*iCH3O2*iNO+ick[,,49]*iHOCH2CH2O2*iNO+
    ick[,,50]*iHOCH2CH2O2*iNO+ick[,,51]*iRN9O2*iNO+ick[,,52]*iCH3CO3*iNO+ick[,,53]*iHOCH2CO3*iNO+ick[,,54]*iRU14O2*iNO+ick[,,55]*iRU14O2*iNO+
    ick[,,56]*iRU12O2*iNO+ick[,,57]*iRU12O2*iNO+ick[,,58]*iRU10O2*iNO+ick[,,59]*iRU10O2*iNO+ick[,,60]*iRU10O2*iNO+ick[,,101]*iCH3NO3+
    ick[,,102]*iOH*iCH3NO3+ick[,,103]*iOH*iHOC2H4NO3+ick[,,104]*iOH*iRN9NO3+ick[,,105]*iOH*iRU14NO3+ick[,,125]*iPAN+
    ick[,,127]*iPHAN+ick[,,128]*iOH*iPAN+ick[,,129]*iOH*iPHAN+ick[,,131]*iRU12PAN+ick[,,133]*iMPAN+ick[,,134]*iOH*iMPAN+ick[,,135]*iOH*iRU12PAN;
  # solver species no.4: O3
  iPoO3 = ick[,,23]*iNO2+ick[,,25]*iNO3;
  # solver species no.5: H2
  iPoH2 = ick[,,39]*iHCHO;
  # solver species no.6: CO
  iPoCO = ick[,,31]*iO3*iC2H4+ick[,,33]*iO3*iC3H6+ick[,,36]*iO3*iC5H8+ick[,,38]*iHCHO+ick[,,39]*iHCHO+ick[,,40]*iCH3CHO+
    ick[,,41]*iOH*iHCHO+ick[,,57]*iRU12O2*iNO+ick[,,88]*iHOCH2CHO+ick[,,90]*iCARB6+ick[,,91]*iUCARB12+
    ick[,,94]*iO3*iUCARB10+ick[,,97]*iOH*iCARB6+ick[,,99]*iO3*iUCARB12+ick[,,128]*iOH*iPAN+ick[,,129]*iOH*iPHAN+ick[,,134]*iOH*iMPAN;
  # solver species no.7: H2O2
  iPoH2O2 = ick[,,11]*iHO2*iHO2+ick[,,12]*iHO2*iHO2+ick[,,95]*iO3*iUCARB10+ick[,,100]*iO3*iUCARB12;
  # solver species no.8: HONO
  iPoHONO = ick[,,13]*iOH*iNO;
  # solver species no.9: HNO3
  iPoHNO3 = ick[,,14]*iOH*iNO2;
  # solver species no.10: HO2NO2
  iPoHO2NO2 = ick[,,17]*iHO2*iNO2;
  # set these units and zero value for the species with no production terms
  # solver species no.11: CH4
  iPoCH4 = array(0.0,c(nbv,nbh));
  # solver species no.12: C2H4
  iPoC2H4 = array(0.0,c(nbv,nbh));
  # solver species no.13: C3H6
  iPoC3H6 = array(0.0,c(nbv,nbh));
  # solver species no.14: C5H8
  iPoC5H8 = array(0.0,c(nbv,nbh));
  # solver species no.15: C2H5OH
  iPoC2H5OH = array(0.0,c(nbv,nbh));
  # The end
  # solver species no.16: CH3NO3
  iPoCH3NO3 = ick[,,61]*iCH3O2*iNO;
  # solver species no.20: CH3OH
  iPoCH3OH = ick[,,75]*iCH3O2;
  # solver species no.22: CH3CO2H
  iPoCH3CO2H = ick[,,34]*iO3*iC3H6;
  # solver species no.23: CH3CHO
  iPoCH3CHO = ick[,,44]*iOH*iC2H5OH+ick[,,51]*iRN9O2*iNO+ick[,,77]*iRN9O2+ick[,,114]*iRN9OOH;
  # solver species no.25: UCARB10
  iPoUCARB10 = ick[,,36]*iO3*iC5H8+ick[,,37]*iO3*iC5H8+ick[,,55]*iRU14O2*iNO+ick[,,81]*iRU14O2+ick[,,110]*iRU14OOH+ick[,,135]*iOH*iRU12PAN;
  # solver species no.26: HCHO
  iPoHCHO = ick[,,31]*iO3*iC2H4+ick[,,32]*iO3*iC2H4+ick[,,33]*iO3*iC3H6+ick[,,34]*iO3*iC3H6+ick[,,43]*iOH*iCH3OH+ick[,,48]*iCH3O2*iNO+2*ick[,,49]*iHOCH2CH2O2*iNO+
    ick[,,51]*iRN9O2*iNO+ick[,,53]*iHOCH2CO3*iNO+ick[,,55]*iRU14O2*iNO+ick[,,59]*iRU10O2*iNO+ick[,,60]*iRU10O2*iNO+ick[,,73]*iCH3O2+ick[,,74]*iCH3O2+ick[,,77]*iRN9O2+ick[,,79]*iHOCH2CO3+ick[,,81]*iRU14O2+ick[,,85]*iRU10O2+ick[,,86]*iRU10O2+
    ick[,,87]*iCARB7+ick[,,88]*iHOCH2CHO+ick[,,89]*iUCARB10+ick[,,94]*iO3*iUCARB10+ick[,,95]*iO3*iUCARB10+ick[,,101]*iCH3NO3+ick[,,102]*iOH*iCH3NO3+ick[,,106]*iCH3OOH+ick[,,108]*iHOCH2CO3H+ick[,,110]*iRU14OOH+
    2*ick[,,113]*iHOC2H4OOH+ick[,,114]*iRN9OOH+ick[,,116]*iOH*iCH3OOH+ick[,,128]*iOH*iPAN+ick[,,129]*iOH*iPHAN;
  # solver species no.28: HOCH2CHO
  iPoHOCH2CHO = ick[,,50]*iHOCH2CH2O2*iNO+ick[,,56]*iRU12O2*iNO+ick[,,58]*iRU10O2*iNO+ick[,,76]*iHOCH2CH2O2+ick[,,82]*iRU12O2+ick[,,83]*iRU12O2+ick[,,84]*iRU10O2+ick[,,91]*iUCARB12+
    ick[,,99]*iO3*iUCARB12+ick[,,100]*iO3*iUCARB12+ick[,,103]*iOH*iHOC2H4NO3+ick[,,111]*iRU12OOH+ick[,,112]*iRU10OOH+ick[,,122]*iOH*iHOC2H4OOH;
  # solver species no.30: HCOOH
  iPoHCOOH = ick[,,32]*iO3*iC2H4+ick[,,37]*iO3*iC5H8;
  # solver species no.31: CARB6
  iPoCARB6 = ick[,,59]*iRU10O2*iNO+ick[,,85]*iRU10O2+ick[,,92]*iOH*iCARB7+ick[,,95]*iO3*iUCARB10+ick[,,100]*iO3*iUCARB12+ick[,,111]*iRU12OOH;
  # solver species no.32: UCARB12
  iPoUCARB12 = ick[,,54]*iRU14O2*iNO+ick[,,80]*iRU14O2+ick[,,105]*iOH*iRU14NO3+ick[,,109]*iRU14OOH+ick[,,119]*iOH*iRU14OOH;
  # solver species no.34: CARB7
  iPoCARB7 = ick[,,57]*iRU12O2*iNO+ick[,,60]*iRU10O2*iNO+ick[,,83]*iRU12O2+ick[,,86]*iRU10O2+ick[,,104]*iOH*iRN9NO3+ick[,,123]*iOH*iRN9OOH+ick[,,134]*iOH*iMPAN;
  # solver species no.36: HOC2H4NO3
  iPoHOC2H4NO3 = ick[,,62]*iHOCH2CH2O2*iNO;
  # solver species no.37: RN9NO3
  iPoRN9NO3 = ick[,,63]*iRN9O2*iNO;
  # solver species no.38: RU14NO3
  iPoRU14NO3 = ick[,,64]*iRU14O2*iNO;
  # solver species no.39: CH3OOH
  iPoCH3OOH = ick[,,65]*iCH3O2*iHO2;
  # solver species no.40: HOC2H4OOH
  iPoHOC2H4OOH = ick[,,66]*iHOCH2CH2O2*iHO2;
  # solver species no.41: RN9OOH
  iPoRN9OOH = ick[,,67]*iRN9O2*iHO2;
  # solver species no.42: CH3CO3H
  iPoCH3CO3H = ick[,,68]*iCH3CO3*iHO2;
  # solver species no.43: HOCH2CO3H
  iPoHOCH2CO3H = ick[,,69]*iHOCH2CO3*iHO2;
  # solver species no.44: RU14OOH
  iPoRU14OOH = ick[,,70]*iRU14O2*iHO2;
  # solver species no.45: RU12OOH
  iPoRU12OOH = ick[,,71]*iRU12O2*iHO2;
  # solver species no.46: RU10OOH
  iPoRU10OOH = ick[,,72]*iRU10O2*iHO2;
  # solver species no.48: PAN
  iPoPAN = ick[,,124]*iCH3CO3*iNO2;
  # solver species no.49: PHAN
  iPoPHAN = ick[,,126]*iHOCH2CO3*iNO2;
  # solver species no.50: RU12PAN
  iPoRU12PAN = ick[,,130]*iRU12O2*iNO2;
  # solver species no.51: MPAN
  iPoMPAN = ick[,,132]*iRU10O2*iNO2;
  
  ## =======================================================================
  ##   *********************** slow reactions **************************  ##
  ##  updated by YQ- on 23/07/2020
  ## =======================================================================
  
  ## Calculation of the chemical production and loss terms
  
  # Solver Species no.1: NO
  iNO = iNO+(iPoNO-iLoNO*iNO)*idtlong;
  
  # Solver Species no.2: NO2
  iNO2 = iNO2+(iPoNO2-iLoNO2*iNO2)*idtlong;
  
  # Solver Species no.4: O3
  iO3 = iO3+(iPoO3-iLoO3*iO3)*idtlong;
  
  # Solver Species no.5: H2
  iH2 = iH2+(iPoH2-iLoH2*iH2)*idtlong;
  
  # Solver Species no.6: CO
  iCO = iCO+(iPoCO-iLoCO*iCO)*idtlong;
  
  # Solver Species no.7: H2O2
  iH2O2 = iH2O2+(iPoH2O2-iLoH2O2*iH2O2)*idtlong;
  
  # Solver Species no.8: HONO
  iHONO = iHONO+(iPoHONO-iLoHONO*iHONO)*idtlong;
  
  # Solver Species no.9: HNO3
  iHNO3 = iHNO3+(iPoHNO3-iLoHNO3*iHNO3)*idtlong;
  
  # Solver Species no.10: HO2NO2
  iHO2NO2 = iHO2NO2+(iPoHO2NO2-iLoHO2NO2*iHO2NO2)*idtlong;
  
  # Solver Species no.11: CH4
  iCH4 = iCH4+(iPoCH4-iLoCH4*iCH4)*idtlong;
  
  # Solver Species no.12: C2H4
  iC2H4 = iC2H4+(iPoC2H4-iLoC2H4*iC2H4)*idtlong;
  
  # Solver Species no.13: C3H6
  iC3H6 = iC3H6+(iPoC3H6-iLoC3H6*iC3H6)*idtlong;
  
  # Solver Species no.14: C5H8
  iC5H8 = iC5H8+(iPoC5H8-iLoC5H8*iC5H8)*idtlong;
  
  # Solver Species no.15: C2H5OH
  iC2H5OH = iC2H5OH+(iPoC2H5OH-iLoC2H5OH*iC2H5OH)*idtlong;
  
  # Solver Species no.16: CH3NO3
  iCH3NO3 = iCH3NO3+(iPoCH3NO3-iLoCH3NO3*iCH3NO3)*idtlong;
  
  # Solver Species no.20: CH3OH
  iCH3OH = iCH3OH+(iPoCH3OH-iLoCH3OH*iCH3OH)*idtlong;
  
  # Solver Species no.22: CH3CO2H
  iCH3CO2H = iCH3CO2H+(iPoCH3CO2H-iLoCH3CO2H*iCH3CO2H)*idtlong;
  
  # Solver Species no.23: CH3CHO
  iCH3CHO = iCH3CHO+(iPoCH3CHO-iLoCH3CHO*iCH3CHO)*idtlong;
  
  # Solver Species no.25: UCARB10
  iUCARB10 = iUCARB10+(iPoUCARB10-iLoUCARB10*iUCARB10)*idtlong;
  
  # Solver Species no.26: HCHO
  iHCHO = iHCHO+(iPoHCHO-iLoHCHO*iHCHO)*idtlong;
  
  # Solver Species no.28: HOCH2CHO
  iHOCH2CHO = iHOCH2CHO+(iPoHOCH2CHO-iLoHOCH2CHO*iHOCH2CHO)*idtlong;
  
  # Solver Species no.30: HCOOH
  iHCOOH = iHCOOH+(iPoHCOOH-iLoHCOOH*iHCOOH)*idtlong;
  
  # Solver Species no.31: CARB6
  iCARB6 = iCARB6+(iPoCARB6-iLoCARB6*iCARB6)*idtlong;
  
  # Solver Species no.32: UCARB12
  iUCARB12 = iUCARB12+(iPoUCARB12-iLoUCARB12*iUCARB12)*idtlong;
  
  # Solver Species no.34: CARB7
  iCARB7 = iCARB7+(iPoCARB7-iLoCARB7*iCARB7)*idtlong;
  
  # Solver Species no.36: HOC2H4NO3
  iHOC2H4NO3 = iHOC2H4NO3+(iPoHOC2H4NO3-iLoHOC2H4NO3*iHOC2H4NO3)*idtlong;
  
  # Solver Species no.37: RN9NO3
  iRN9NO3 = iRN9NO3+(iPoRN9NO3-iLoRN9NO3*iRN9NO3)*idtlong;
  
  # Solver Species no.38: RU14NO3
  iRU14NO3 = iRU14NO3+(iPoRU14NO3-iLoRU14NO3*iRU14NO3)*idtlong;
  
  # Solver Species no.39: CH3OOH
  iCH3OOH = iCH3OOH+(iPoCH3OOH-iLoCH3OOH*iCH3OOH)*idtlong;
  
  # Solver Species no.40: HOC2H4OOH
  iHOC2H4OOH = iHOC2H4OOH+(iPoHOC2H4OOH-iLoHOC2H4OOH*iHOC2H4OOH)*idtlong;
  
  # Solver Species no.41: RN9OOH
  iRN9OOH = iRN9OOH+(iPoRN9OOH-iLoRN9OOH*iRN9OOH)*idtlong;
  
  # Solver Species no.42: CH3CO3H
  iCH3CO3H = iCH3CO3H+(iPoCH3CO3H-iLoCH3CO3H*iCH3CO3H)*idtlong;
  
  # Solver Species no.43: HOCH2CO3H
  iHOCH2CO3H = iHOCH2CO3H+(iPoHOCH2CO3H-iLoHOCH2CO3H*iHOCH2CO3H)*idtlong;
  
  # Solver Species no.44: RU14OOH
  iRU14OOH = iRU14OOH+(iPoRU14OOH-iLoRU14OOH*iRU14OOH)*idtlong;
  
  # Solver Species no.45: RU12OOH
  iRU12OOH = iRU12OOH+(iPoRU12OOH-iLoRU12OOH*iRU12OOH)*idtlong;
  
  # Solver Species no.46: RU10OOH
  iRU10OOH = iRU10OOH+(iPoRU10OOH-iLoRU10OOH*iRU10OOH)*idtlong;
  
  # Solver Species no.48: PAN
  iPAN = iPAN+(iPoPAN-iLoPAN*iPAN)*idtlong;
  
  # Solver Species no.49: PHAN
  iPHAN = iPHAN+(iPoPHAN-iLoPHAN*iPHAN)*idtlong;
  
  # Solver Species no.50: RU12PAN
  iRU12PAN = iRU12PAN+(iPoRU12PAN-iLoRU12PAN*iRU12PAN)*idtlong;
  
  # Solver Species no.51: MPAN
  iMPAN = iMPAN+(iPoMPAN-iLoMPAN*iMPAN)*idtlong
  
  ## =======================================================================
  ##   *********************** fast reactions **************************  ##
  ##  updated by YQ- on 23/07/2020
  ## ======================================================================= 
  while (counter < nshort) 
  {
    ##  ******************* calculation of loss term ******************** ##
    # Solver Species no.3: NO3     
    iLoNO3 = ick[,,4]*iNO+ick[,,15]*iOH+ick[,,24]+ick[,,25];
    # Solver Species no.17: OH
    iLoOH = ick[,,5]*iO3+ick[,,6]*iH2+ick[,,7]*iCO+ick[,,8]*iH2O2+ick[,,10]*iHO2+ick[,,13]*iNO+ick[,,14]*iNO2+ick[,,15]*iNO3+
      ick[,,19]*iHO2NO2+ick[,,20]*iHONO+ick[,,21]*iHNO3+ick[,,28]*iCH4+ick[,,29]*iC2H4+ick[,,30]*iC3H6+ick[,,35]*iC5H8+
      ick[,,41]*iHCHO+ick[,,42]*iCH3CHO+ick[,,43]*iCH3OH+ick[,,44]*iC2H5OH+ick[,,45]*iC2H5OH+ick[,,46]*iHCOOH+ick[,,47]*iCH3CO2H+ick[,,92]*iCARB7+ick[,,93]*iUCARB10+
      ick[,,96]*iHOCH2CHO+ick[,,97]*iCARB6+ick[,,98]*iUCARB12+ick[,,102]*iCH3NO3+ick[,,103]*iHOC2H4NO3+ick[,,104]*iRN9NO3+ick[,,105]*iRU14NO3+
      ick[,,115]*iCH3OOH+ick[,,116]*iCH3OOH+
      ick[,,117]*iCH3CO3H+ick[,,118]*iHOCH2CO3H+ick[,,119]*iRU14OOH+ick[,,120]*iRU12OOH+ick[,,121]*iRU10OOH+ick[,,122]*iHOC2H4OOH+ick[,,123]*iRN9OOH+
      ick[,,128]*iPAN+ick[,,129]*iPHAN+ick[,,134]*iMPAN+ick[,,135]*iRU12PAN;
    # Solver Species no.18: HO2
    iLoHO2 = ick[,,9]*iO3+ick[,,10]*iOH+2*ick[,,11]*iHO2+2*ick[,,12]*iHO2+ick[,,16]*iNO+
      ick[,,17]*iNO2+
      ick[,,65]*iCH3O2+ick[,,66]*iHOCH2CH2O2+ick[,,67]*iRN9O2+
      ick[,,68]*iCH3CO3+ick[,,69]*iHOCH2CO3+ick[,,70]*iRU14O2+ick[,,71]*iRU12O2+ick[,,72]*iRU10O2;
    # Solver Species no.19: CH3O2
    iLoCH3O2 = ick[,,48]*iNO+ick[,,61]*iNO+ick[,,65]*iHO2+ick[,,73]+ick[,,74]+ick[,,75];
    # Solver Species no.21: HOCH2CH2O2
    iLoHOCH2CH2O2 = ick[,,49]*iNO+ick[,,50]*iNO+ick[,,62]*iNO+ick[,,66]*iHO2+ick[,,76];
    # Solver Species no.24: RU14O2
    iLoRU14O2 = ick[,,54]*iNO+ick[,,55]*iNO+ick[,,64]*iNO+ick[,,70]*iHO2+ick[,,80]+ick[,,81];
    # Solver Species no.27: CH3CO3
    iLoCH3CO3 = ick[,,52]*iNO+ick[,,68]*iHO2+ick[,,78]+ick[,,124]*iNO2;
    # Solver Species no.29: RN9O2
    iLoRN9O2 = ick[,,51]*iNO+ick[,,63]*iNO+ick[,,67]*iHO2+ick[,,77];
    # Solver Species no.33: RU12O2
    iLoRU12O2 = ick[,,56]*iNO+ick[,,57]*iNO+ick[,,71]*iHO2+ick[,,82]+ick[,,83]+ick[,,130]*iNO2;
    # Solver Species no.35: RU10O2
    iLoRU10O2 = ick[,,58]*iNO+ick[,,59]*iNO+ick[,,60]*iNO+ick[,,72]*iHO2+ick[,,84]+ick[,,85]+ick[,,86]+ick[,,132]*iNO2;
    # Solver Species no.47: HOCH2CO3
    iLoHOCH2CO3 = ick[,,53]*iNO+ick[,,69]*iHO2+ick[,,79]+ick[,,126]*iNO2;
    
    ##  ****************** calculation of production term ******************* ##
    # solver species no.3: NO3
    iPoNO3=ick[,,21]*iOH*iHNO3+ick[,,136]*iNO2*iO3;
    # solver species no.17: OH
    iPoOH = 2*ick[,,1]*iO3+ick[,,9]*iHO2*iO3+ick[,,16]*iHO2*iNO+
      2*ick[,,22]*iH2O2+ick[,,26]*iHONO+ick[,,27]*iHNO3+ick[,,31]*iO3*iC2H4+ick[,,33]*iO3*iC3H6+
      ick[,,36]*iO3*iC5H8+ick[,,94]*iO3*iUCARB10+
      ick[,,99]*iO3*iUCARB12+ick[,,106]*iCH3OOH+ick[,,107]*iCH3CO3H+
      ick[,,108]*iHOCH2CO3H+ick[,,109]*iRU14OOH+ick[,,110]*iRU14OOH+ick[,,111]*iRU12OOH+ick[,,112]*iRU10OOH+ick[,,113]*iHOC2H4OOH+ick[,,114]*iRN9OOH+ick[,,116]*iOH*iCH3OOH+
      ick[,,119]*iOH*iRU14OOH+ick[,,122]*iOH*iHOC2H4OOH+ick[,,123]*iOH*iRN9OOH;
    
    # solver species no.18: HO2
    iPoHO2 = ick[,,5]*iOH*iO3+ick[,,6]*iOH*iH2+ick[,,7]*iOH*iCO+ick[,,8]*iOH*iH2O2+ick[,,15]*iOH*iNO3+
      ick[,,18]*iHO2NO2+ick[,,31]*iO3*iC2H4+ick[,,36]*iO3*iC5H8+2*ick[,,38]*iHCHO+ick[,,40]*iCH3CHO+ick[,,41]*iOH*iHCHO+ick[,,43]*iOH*iCH3OH+ick[,,44]*iOH*iC2H5OH+ick[,,46]*iHCOOH*iOH+ick[,,48]*iCH3O2*iNO+
      ick[,,49]*iHOCH2CH2O2*iNO+ick[,,50]*iHOCH2CH2O2*iNO+ick[,,51]*iRN9O2*iNO+ick[,,53]*iHOCH2CO3*iNO+ick[,,54]*iRU14O2*iNO+ick[,,55]*iRU14O2*iNO+ick[,,57]*iRU12O2*iNO+ick[,,59]*iRU10O2*iNO+ick[,,60]*iRU10O2*iNO+
      ick[,,73]*iCH3O2+ick[,,76]*iHOCH2CH2O2+ick[,,77]*iRN9O2+ick[,,79]*iHOCH2CO3+ick[,,80]*iRU14O2+ick[,,81]*iRU14O2+ick[,,83]*iRU12O2+
      ick[,,85]*iRU10O2+ick[,,86]*iRU10O2+ick[,,87]*iCARB7+2*ick[,,88]*iHOCH2CHO+ick[,,89]*iUCARB10+ick[,,90]*iCARB6+ick[,,91]*iUCARB12+ick[,,92]*iOH*iCARB7+ick[,,101]*iCH3NO3+ick[,,106]*iCH3OOH+
      ick[,,108]*iHOCH2CO3H+ick[,,109]*iRU14OOH+ick[,,110]*iRU14OOH+ick[,,111]*iRU12OOH+ick[,,113]*iHOC2H4OOH+ick[,,114]*iRN9OOH;
    
    # solver species no.19: CH3O2
    iPoCH3O2 = ick[,,28]*iOH*iCH4+ick[,,33]*iO3*iC3H6+ick[,,40]*iCH3CHO+ick[,,47]*iCH3CO2H*iOH+ick[,,52]*iCH3CO3*iNO+ick[,,78]*iCH3CO3+ick[,,107]*iCH3CO3H+ick[,,115]*iOH*iCH3OOH;
    #solver species no.21: HOCH2CH2O2
    iPoHOCH2CH2O2 = ick[,,29]*iOH*iC2H4+ick[,,45]*iOH*iC2H5OH;
    # solver species no.24: RU14O2
    iPoRU14O2 = ick[,,35]*iOH*iC5H8;
    # solver species no.27: CH3CO3
    iPoCH3CO3 = ick[,,42]*iOH*iCH3CHO+ick[,,56]*iRU12O2*iNO+ick[,,58]*iRU10O2*iNO+ick[,,82]*iRU12O2+ick[,,84]*iRU10O2+ick[,,87]*iCARB7+ick[,,89]*iUCARB10+
      ick[,,90]*iCARB6+ick[,,91]*iUCARB12+ick[,,94]*iO3*iUCARB10+ick[,,97]*iOH*iCARB6+ick[,,99]*iO3*iUCARB12+ick[,,112]*iRU10OOH+ick[,,117]*iOH*iCH3CO3H+ick[,,125]*iPAN;
    # solver species no.29: RN9O2
    iPoRN9O2 = ick[,,30]*iOH*iC3H6;
    # solver species no.33: RU12O2
    iPoRU12O2 = ick[,,98]*iOH*iUCARB12+ick[,,120]*iOH*iRU12OOH+ick[,,131]*iRU12PAN;
    # solver species no.35: RU10O2
    iPoRU10O2 = ick[,,93]*iOH*iUCARB10+ick[,,121]*iOH*iRU10OOH+ick[,,133]*iMPAN;
    # solver species no.47: HOCH2CO3
    iPoHOCH2CO3 = ick[,,96]*iOH*iHOCH2CHO+ick[,,118]*iOH*iHOCH2CO3H+ick[,,127]*iPHAN;
    
    ##  ****************** calculation of fast chemsitry ******************* ##
    # Solver Species no.3: NO3
    iNO3 = iPoNO3/iLoNO3+(iNO3-iPoNO3/iLoNO3)/(1+dtchem*iLoNO3+0.5*dtchem*dtchem*iLoNO3*iLoNO3);
    # Solver Species no.17: OH
    iOH = iPoOH/iLoOH+(iOH-iPoOH/iLoOH)/(1+dtchem*iLoOH+0.5*dtchem*dtchem*iLoOH*iLoOH);
    # Solver Species no.18: HO2
    iHO2 = iPoHO2/iLoHO2+(iHO2-iPoHO2/iLoHO2)/(1+dtchem*iLoHO2+0.5*dtchem*dtchem*iLoHO2*iLoHO2);
    # Solver Species no.19: CH3O2
    iCH3O2 = iPoCH3O2/iLoCH3O2+(iCH3O2-iPoCH3O2/iLoCH3O2)/(1+dtchem*iLoCH3O2+0.5*dtchem*dtchem*iLoCH3O2*iLoCH3O2);
    # Solver Species no.21: HOCH2CH2O2
    iHOCH2CH2O2 = iPoHOCH2CH2O2/iLoHOCH2CH2O2+(iHOCH2CH2O2-iPoHOCH2CH2O2/iLoHOCH2CH2O2)/(1+dtchem*iLoHOCH2CH2O2+0.5*dtchem*dtchem*iLoHOCH2CH2O2*iLoHOCH2CH2O2);
    # Solver Species no.24: RU14O2
    iRU14O2 = iPoRU14O2/iLoRU14O2+(iRU14O2-iPoRU14O2/iLoRU14O2)/(1+dtchem*iLoRU14O2+0.5*dtchem*dtchem*iLoRU14O2*iLoRU14O2);
    # Solver Species no.27: CH3CO3
    iCH3CO3 = iPoCH3CO3/iLoCH3CO3+(iCH3CO3-iPoCH3CO3/iLoCH3CO3)/(1+dtchem*iLoCH3CO3+0.5*dtchem*dtchem*iLoCH3CO3*iLoCH3CO3);
    # Solver Species no.29: RN9O2
    iRN9O2 = iPoRN9O2/iLoRN9O2+(iRN9O2-iPoRN9O2/iLoRN9O2)/(1+dtchem*iLoRN9O2+0.5*dtchem*dtchem*iLoRN9O2*iLoRN9O2);
    # Solver Species no.33: RU12O2
    iRU12O2 = iPoRU12O2/iLoRU12O2+(iRU12O2-iPoRU12O2/iLoRU12O2)/(1+dtchem*iLoRU12O2+0.5*dtchem*dtchem*iLoRU12O2*iLoRU12O2);
    # Solver Species no.35: RU10O2
    iRU10O2 = iPoRU10O2/iLoRU10O2+(iRU10O2-iPoRU10O2/iLoRU10O2)/(1+dtchem*iLoRU10O2+0.5*dtchem*dtchem*iLoRU10O2*iLoRU10O2);
    # Solver Species no.47: HOCH2CO3
    iHOCH2CO3 = iPoHOCH2CO3/iLoHOCH2CO3+(iHOCH2CO3-iPoHOCH2CO3/iLoHOCH2CO3)/(1+dtchem*iLoHOCH2CO3+0.5*dtchem*dtchem*iLoHOCH2CO3*iLoHOCH2CO3);
    #
    counter = counter + 1
  }
  
  ## ********************** chemical outputs  *************************** ##
  #  Initial species NO 01-10
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
  #  Initial species NO 11-20
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
  #  Initial species NO 21-30
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
  #  Initial species NO 31-40
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
  #  Initial species NO 41-51
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
  return(Carr)
}
## **************************  COMPLETE!  ***************************** ##
