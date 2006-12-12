//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// author: Vladimir.Grichine@cern.ch
//
// Implements data from: Barashenkov V.S., Nucleon-Nucleus Cross Section,
// Preprint JINR P2-89-770, p. 12, Dubna 1989 (scanned version from KEK)
// Based on G. Folger version of G4PiNuclearCrossSection class
//

#include "G4NucleonNuclearCrossSection.hh"
#include "G4HadronicException.hh"
#include "G4HadTmpUtil.hh"
#include "G4ping.hh"

// Group 1: He, Be, C for 44 energies  

const G4double G4NucleonNuclearCrossSection::e1[44] =     
{
  0.014, 0.015, 0.017, .02, 0.022, 0.025, 0.027, 0.03, 0.035, .04, 0.045, 0.05, .06, 0.07, 
  .08, 0.09,  .1, .12, .14, .15, .16, .18, .20, .25, .30, .35, .4 , 0.5, 0.6, 0.7,  0.8,  
  0.9,   1, 1.5,   2,   3,   5,  7, 10,   
  20,   50,  100,  500, 1000
};
const G4double G4NucleonNuclearCrossSection::he_m_t[44] =   
{ 
  1090, 1020, 915, 800, 710, 640, 600, 560, 500, 440, 390, 360, 295, 256, 220, 192, 
  168, 136, 120, 116, 114, 110, 107, 104, 106, 108, 110, 120, 126, 135, 140, 144, 146, 
  148, 152, 150, 146, 142, 138, 132, 129, 126, 127, 128  
};
const G4double G4NucleonNuclearCrossSection::he_m_in[44] =  
{ 
  0, 5, 10, 20, 35, 55, 70, 80, 90, 105, 115, 115, 100, 90, 86, 84, 84, 82, 80, 80, 80, 80, 
  79, 78, 80, 84, 88, 94, 100, 105, 108, 108, 108, 112, 114, 114, 112, 110, 108, 106, 104, 
  101, 102, 102
};
const G4double G4NucleonNuclearCrossSection::he_t_in[44] =  
{ 
  0, 2, 3, 13, 30, 50, 65, 77, 90, 105, 115, 115, 100, 90, 86, 84, 84, 82, 80, 80, 80, 80, 
  79, 78, 80, 84, 88, 94, 100, 105, 108, 108, 108, 112, 114, 114, 112, 110, 108, 106, 104, 
  101, 102, 102
};
const G4double G4NucleonNuclearCrossSection::be_m_t[44] = 
{
  1490, 1460, 1400, 1350, 1270, 1200, 1160, 1100, 1000, 910, 810, 740, 625, 575, 455, 406, 
  365, 310, 275, 262, 255, 240, 235, 225, 225, 230, 238, 252, 270, 282, 288, 290, 294, 303, 
  303, 300, 292, 284, 277, 267, 263, 264, 268, 268 
};
const G4double G4NucleonNuclearCrossSection::be_m_in[44] =
{ 
  650, 640, 617, 595, 555, 520, 495, 470, 430, 385, 350, 320, 270, 250, 210, 190, 185, 178, 
  175, 175, 175, 175, 175, 170, 170, 172, 176, 184, 194, 200, 209, 213, 214, 216, 216, 212, 
  210, 210, 210, 210, 210, 210, 210, 210
};
const G4double G4NucleonNuclearCrossSection::be_p_in[44] =
{ 
  490, 540, 580, 545, 525, 495, 470, 450, 420, 370, 340, 310, 262, 242, 205, 185, 180, 175, 
  172, 175, 175, 175, 175, 170, 170, 172, 176, 184, 194, 200, 209, 213, 214, 216, 216, 212, 
  210, 210, 210, 210, 210, 210, 210, 210
};
const G4double G4NucleonNuclearCrossSection::c_m_t[44] = 
{
  1240, 1370, 1450, 1455, 1445, 1385, 1345, 1290, 1210, 1110, 1020, 940, 800, 700, 604, 530, 
  475, 396, 350, 336, 320, 303, 294, 280, 280, 286, 296, 314, 330, 344, 356, 360, 364, 384, 
  388, 384, 364, 352, 344, 330, 324, 324, 332, 332
};
const G4double G4NucleonNuclearCrossSection::c_m_in[44] =
{
  590, 570, 542, 510, 500, 460, 445, 430, 395, 380, 350, 330, 295, 270, 255, 240, 228, 222, 
  216, 216, 210, 210, 210, 208, 210, 214, 216, 228, 240, 248, 254, 257, 260, 262, 260, 256, 
  252, 252, 250, 250, 248, 248, 248, 248
};
const G4double G4NucleonNuclearCrossSection::c_p_in[44] =
{ 
  310, 330, 400, 440, 450, 435, 430, 420, 385, 370, 340, 320, 288, 263, 249, 234, 222, 216, 
  210, 211, 205, 208, 210, 208, 210, 214, 216, 228, 240, 248, 254, 257, 260, 262, 260, 256, 
  252, 252, 250, 250, 248, 248, 248, 248
};

// Group 2: N, O, Na for 44 energies (e1=e2)

const G4double G4NucleonNuclearCrossSection::e2[44] =    
{
 0.014, 0.015, 0.017, .02, 0.022, 0.025, 0.027, 0.03, 0.035, .04, 0.045, 0.05, .06, 0.07, 
  .08, 0.09,  .1, .12, .14, .15, .16, .18, .20, .25, .30, .35, .4 , 0.5, 0.6, 0.7,  0.8,  
  0.9,   1, 1.5,   2,   3,   5,  7, 10,   
  20,   50,  100,  500, 1000  
};
const G4double G4NucleonNuclearCrossSection::n_m_t[44] = 
{
  1420,1480, 1537, 1550, 1525, 1500, 1480, 1425, 1340, 1260, 1175, 1090, 930, 805, 690, 612, 
  552, 462, 402, 384, 372, 350, 345, 326, 324, 328, 336, 356, 372, 388, 400, 408, 415, 430, 
  435, 432, 415, 402, 390, 375, 367, 370, 382, 385
};
const G4double G4NucleonNuclearCrossSection::n_m_in[44] =
{
  680, 665, 625, 580, 562, 525, 510, 485, 450, 435, 410, 387, 340, 310, 290, 280, 276, 274, 
  260, 258, 254, 247, 245, 240, 240, 244, 250, 260, 268, 275, 280, 285, 290, 295, 300, 294, 
  292, 290, 285, 285, 282, 282, 282, 282
};
const G4double G4NucleonNuclearCrossSection::n_p_in[44] =
{ 
  420, 440, 470, 490, 497, 500, 480, 462, 440, 425, 400, 377, 333, 303, 284, 274, 270, 268, 
  254, 252, 247, 245, 245, 240, 240, 244, 250, 260, 268, 275, 280, 285, 290, 295, 300, 294, 
  292, 290, 285, 285, 282, 282, 282, 282
};
const G4double G4NucleonNuclearCrossSection::e3[31] =     
{
  .02, .04, .06, .08,  0.1,  .12,  .14,  .16, 0.18, 0.2, .22, .25, 0.3, .35, 0.4, 
.45, 0.5, 0.6, 0.7, 0.8, 0.9,   1,   2,   3,   5,  10,  20,  50, 100, 500, 1000
};
const G4double G4NucleonNuclearCrossSection::o_m_t[31] =  
{
  280, 360, 500, 685,  812,  861,  870,  865,  835, 800, 755, 700, 600, 537, 493, 
  468, 441, 436, 443, 449, 460, 463, 432, 385, 350, 325, 312, 307, 303, 303,  303
};
const G4double G4NucleonNuclearCrossSection::o_m_in[31] = 
{
  190, 207, 300, 420,  500,  540,  550,  542,  520, 490, 460, 423, 360, 339, 321, 
314, 312, 314, 319, 324, 328, 330, 300, 275, 250, 240, 229, 225, 222, 222,  222
};
const G4double G4NucleonNuclearCrossSection::o_p_t[20] =  
{
  170, 240, 390, 570,  740,  818,  830,  822,  800, 765, 725, 675, 585, 525, 483, 
458, 444, 447, 453, 449
};
const G4double G4NucleonNuclearCrossSection::o_p_in[20] = 
{
  100, 145, 240, 340,  470,  518,  530,  522,  505, 477, 448, 412, 350, 330, 316, 310, 
308, 311, 317, 324
};
const G4double G4NucleonNuclearCrossSection::na_m_t[31] = 
{
  450, 545, 705, 910, 1020, 1075, 1087, 1080, 1042, 987, 943, 885, 790, 700, 650, 
610, 585, 575, 585, 595, 600, 610, 556, 524, 494, 458, 445, 429, 427, 427,  427
};
const G4double G4NucleonNuclearCrossSection::na_m_in[31] =
{
  275, 315, 413, 545,  620,  660,  670,  662,  630, 593, 570, 520, 465, 420, 410, 
395, 390, 400, 410, 418, 420, 422, 372, 348, 330, 320, 310, 294, 292, 292,  292
};
const G4double G4NucleonNuclearCrossSection::na_p_t[22] = 
{
  210, 320, 530, 795,  960, 1035, 1050, 1040, 1007, 957, 918, 865, 773, 685, 636, 
  598, 575, 565, 578, 590, 598, 610
};
const G4double G4NucleonNuclearCrossSection::na_p_in[22] =
{
  115, 210, 340, 495,  585,  630,  645,  637,  605, 572, 550, 505, 455, 410, 401, 388, 
  383, 393, 405, 414, 418, 422
};
const G4double G4NucleonNuclearCrossSection::e3_1[31] =  
{ 
  0.02, .04,  .06,  .08,   .1,  .12,  .14,  .16,  .18,   .2,  .22,  .25,   .3,  
  .35,   .4, .45,  .5,  .6,  .7,  .8,  .9,  1.,  2.,  3.,  5., 10., 20., 50., 
100.,500.,1000.
};
const G4double G4NucleonNuclearCrossSection::al_m_t[31] = 
{ 
  532, 637,  832, 1057, 1207, 1230, 1210, 1174, 1133, 1095, 1038,  970,  890,  807,  
  750, 710, 675, 665, 670, 673, 678, 682, 618, 574, 546, 520, 507, 495,  488, 488,  488
};
const G4double G4NucleonNuclearCrossSection::al_m_in[31] = 
{
  300, 360,  495,  665,  750,  765,  750,  730,  700,  660,  615,  570,  520,  490,  
  470, 450, 448, 450, 450, 452, 456, 460, 408, 392, 376, 356, 347, 338,  332, 332,  332
};
const G4double G4NucleonNuclearCrossSection::al_p_t[21] = 
{ 
  225, 350,  616,  945, 1122, 1175, 1157, 1128, 1088, 1045,  988,  935,  870,  787,  
  730, 690, 660, 652, 660, 668, 678
};
const G4double G4NucleonNuclearCrossSection::al_p_in[21] = 
{
  120, 238,  390,  610,  712,  735,  720,  703,  655,  635,  590,  550,  505,  475,  
455, 438, 440, 445, 445, 450, 456
};
const G4double G4NucleonNuclearCrossSection::ca_m_t[31] = 
{ 
  800, 980, 1240, 1460, 1570, 1600, 1580, 1535, 1475, 1425, 1375, 1295, 1200, 1083, 
  1000, 948, 915, 895, 900, 908, 915, 922, 856, 795, 740, 705, 682, 660,  660, 660,  660
};
const G4double G4NucleonNuclearCrossSection::ca_m_in[31] = 
{
  470, 550,  620,  860,  955,  980,  960,  920,  860,  820,  780,  740,  665,  637,  
  615, 600, 590, 590, 600, 608, 610, 615, 550, 525, 510, 488, 470, 450,  450, 450,  450
};
const G4double G4NucleonNuclearCrossSection::ca_p_t[23] = 
{ 
  275, 445,  790, 1195, 1440, 1485, 1475, 1435, 1385, 1335, 1295, 1245, 1160, 1050,  
  970, 923, 895, 877, 887, 897, 904, 913, 855
};
const G4double G4NucleonNuclearCrossSection::ca_p_in[23] = 
{
  160, 315,  500,  745,  870,  905,  900,  860,  810,  770,  740,  710,  640,  617,  
  595, 585, 575, 575, 590, 600, 602, 608, 510
}; // last number is 500 in org, changed to make things smooth.
const G4double G4NucleonNuclearCrossSection::e4[32] =     
{
  0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.14, 0.16, 0.18,  0.2, 0.22, 0.25,  0.3, 0.35,  
  0.4, 0.45,  0.5, 0.55,  0.6,  0.7,  0.8,  0.9,    1,    2,    3,    5,   10,   20,   
  50,  100,  500, 1000
};
const G4double G4NucleonNuclearCrossSection::fe_m_t[32] = 
{
  1175, 1363, 1670, 1950, 2050, 2040, 1975, 1886, 1834, 1773, 1720, 1635, 1474, 1380, 
  1269, 1225, 1182, 1162, 1159, 1162, 1178, 1190, 1197, 1102, 1135,  975,  945,  925,  
  905,  905,  905,  905
};
const G4double G4NucleonNuclearCrossSection::fe_m_in[32] = 
{
  625, 725,   910, 1180, 1275, 1250, 1200, 1150, 1100, 1040,  995,  925,  825,  810,  
  780,  760,  745,  740,  740,  740,  750,  760,  765,  690,  660,  635,  615,  600,  
  585,  585,  585,  585
};
const G4double G4NucleonNuclearCrossSection::fe_p_t[25] =  
{
  330, 575,  1010, 1500, 1837, 1875, 1820, 1751, 1691, 1636, 1690, 1450, 1396, 1305, 
  1219, 1190, 1148, 1138, 1134, 1144, 1163, 1175, 1183, 1198, 1135
};
const G4double G4NucleonNuclearCrossSection::fe_p_in[25] = 
{
  210, 410,   707, 1010, 1125, 1150, 1100, 1070, 1010,  960,  920,  776,  780,  760,  
  750,  740,  720,  725,  725,  730,  740,  750,  755,  690,  660
};
const G4double G4NucleonNuclearCrossSection::cu_m_t[32] = 
{
  1400, 1600, 1875, 2088, 2200, 2220, 2175, 2125, 2075, 2012, 1950, 1855, 1670, 1530, 
  1430, 1370, 1315, 1315, 1315, 1330, 1345, 1360, 1365, 1250, 1185, 1128, 1070, 1035, 
  1010, 1010, 1010, 1010
};
const G4double G4NucleonNuclearCrossSection::cu_m_in[32] = 
{
  725, 840,  1020, 1200, 1295, 1300, 1267, 1240, 1213, 1175, 1125, 1042,  950,  900,  
  860,  840,  830,  832,  835,  840,  850,  860,  865,  785,  735,  705,  680,  650,  
  630,  630,  630,  630
};
const G4double G4NucleonNuclearCrossSection::cu_p_t[25] =  
{
  355, 605,  1120, 1630, 1940, 2010, 2010, 1980, 1925, 1895, 1830, 1730, 1585, 1490, 
  1400, 1340, 1290, 1290, 1290, 1310, 1330, 1345, 1350, 1240, 1185
};
const G4double G4NucleonNuclearCrossSection::cu_p_in[25] = 
{
  230, 425,  780,  1025, 1155, 1190, 1190, 1180, 1125, 1100, 1050, 1000,  900,  870,  
  835,  815,  810,  812,  815,  825,  840,  850,  855,  780,  735
};
const G4double G4NucleonNuclearCrossSection::e5[34] =     
{
  0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.12, 0.14, 0.16, 0.18,  0.2, 0.22, 
  0.25,  0.3, 0.35,  0.4, 0.45,  0.5,  0.6,  0.7,  0.8,  0.9,    1,    2,    3,    5,   
  10,   20,   50,  100,  500, 1000
};
const G4double G4NucleonNuclearCrossSection::mo_m_t[34] = 
{
  2430, 2610, 2710, 2790, 2880, 2940, 2965, 2970, 2970, 2920, 2840, 2720, 2570, 2500, 
  2365, 2200, 2050, 1926, 1825, 1768, 1749, 1750, 1778, 1789, 1808, 1690, 1645, 1530, 
  1492, 1450, 1425, 1425, 1425, 1425
};
const G4double G4NucleonNuclearCrossSection::mo_m_in[34] = 
{
  925, 1125, 1250, 1375, 1500, 1600, 1680, 1750, 1770, 1730, 1660, 1580, 1500, 1450, 
  1330, 1250, 1190, 1140, 1100, 1075, 1075, 1070, 1088, 1095, 1110, 1035, 1005,  940,  
  917,  880,  860,  860,  860,  860
};
const G4double G4NucleonNuclearCrossSection::mo_p_t[27] =  
{
  410,  730, 1110, 1530, 1920, 2200, 2385, 2520, 2600, 2630, 2575, 2470, 2320, 2285, 
  2185, 2053, 1945, 1852, 1776, 1719, 1710, 1716, 1746, 1759, 1778, 1675, 1645
};
const G4double G4NucleonNuclearCrossSection::mo_p_in[27] = 
{
  270,  540,  825,  975, 1140, 1285, 1400, 1480, 1555, 1580, 1525, 1470, 1360, 1340, 
  1255, 1160, 1120, 1085, 1060, 1045, 1045, 1045, 1065, 1075, 1090, 1025, 1005
};
const G4double G4NucleonNuclearCrossSection::cd_m_t[34] = 
{
  3060, 3125, 3170, 3220, 3255, 3280, 3290, 3260, 3270, 3200, 3120, 3080, 3090, 2920, 
  2810, 2640, 2362, 2230, 2115, 2050, 2020, 2025, 2040, 2070, 2100, 1900, 1795, 1740, 
  1675, 1645, 1625, 1620, 1620, 1620
};
const G4double G4NucleonNuclearCrossSection::cd_m_in[34]= 
{
  1025, 1275, 1440, 1625, 1740, 1800, 1880, 1920, 1980, 1920, 1850, 1810, 1720, 1650, 
  1560, 1450, 1330, 1290, 1245, 1210, 1200, 1200, 1205, 1205, 1230, 1130, 1085, 1060, 
  1000,  985,  975,  970,  970,  970
};
const G4double G4NucleonNuclearCrossSection::cd_p_t[28] =  
{
  455,  780, 1170, 1700, 2120, 2400, 2600, 2720, 2820, 2840, 2800, 2760, 2720, 2640, 
  2560, 2450, 2252, 2130, 2035, 1985, 1970, 1975, 2005, 2035, 2070, 1880, 1795, 1740
};
const G4double G4NucleonNuclearCrossSection::cd_p_in[28] = 
{
  310,  580,  880, 1060, 1270, 1400, 1530, 1610, 1660, 1680, 1640, 1600, 1560, 1500, 
  1430, 1330, 1280, 1230, 1200, 1180, 1170, 1175, 1180, 1180, 1210, 1120, 1085, 1060
};
const G4double G4NucleonNuclearCrossSection::e6[35] =      
{
  0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.12, 0.14, 0.16, 0.18,  0.2, 0.22, 
  0.25,  0.3, 0.35,  0.4, 0.45,  0.5, 0.55,  0.6,  0.7,  0.8,  0.9,    1,    2,    3,    
  5,   10,   20,   50,  100,  500, 1000
};
const G4double G4NucleonNuclearCrossSection::sn_m_t[35] =  
{
  3000, 3180, 3250, 3300, 3300, 3410, 3470, 3450, 3410, 3350, 3280, 3200, 3120, 3050, 
  2900, 2630, 2500, 2325, 2190, 2100, 2060, 2055, 2055, 2055, 2067, 2085, 2000, 1900, 
  1835, 1770, 1720, 1700, 1695, 1695, 1695
};
const G4double G4NucleonNuclearCrossSection::sn_m_in[35] = 
{
  1050, 1350, 1520, 1650, 1800, 1980, 2070, 2120, 2090, 2050, 1980, 1920, 1830, 1770, 
  1670, 1500, 1435, 1350, 1300, 1230, 1220, 1235, 1235, 1235, 1237, 1240, 1160, 1120, 
  1090, 1065, 1040, 1020, 1015, 1015, 1015
};
const G4double G4NucleonNuclearCrossSection::sn_p_t[29] =  
{ 
  465,  800, 1200, 1760, 2170, 2480, 2730, 2885, 2970, 2980, 2970, 2890, 2840, 2790, 
  2620, 2450, 2335, 2205, 2080, 2020, 2010, 1990, 1990, 2015, 2030, 2045, 1980, 1890, 1835
};
const G4double G4NucleonNuclearCrossSection::sn_p_in[29] = 
{ 
  315,  590,  880, 1220, 1460, 1580, 1700, 1770, 1810, 1810, 1800, 1730, 1680, 1630, 
  1530, 1400, 1335, 1270, 1210, 1180, 1190, 1190, 1190, 1205, 1210, 1210, 1150, 1115, 1090
};
const G4double G4NucleonNuclearCrossSection::w_m_t[35] =   
{
  5200, 5115, 5025, 4975, 4900, 4850, 4780, 4725, 4600, 4490, 4355, 4255, 4125, 4040, 
  3830, 3580, 3330, 3110, 2955, 2860, 2852, 2845, 2885, 2900, 2915, 2940, 2800, 2660, 
  2570, 2490, 2460, 2425, 2420, 2420, 2420
};
const G4double G4NucleonNuclearCrossSection::w_m_in[35] =  
{
  1450, 1850, 2100, 2350, 2550, 2700, 2825, 2900, 2850, 2750, 2630, 2525, 2400, 2300, 
  2200, 2070, 1880, 1770, 1715, 1680, 1680, 1680, 1685, 1690, 1700, 1720, 1635, 1560, 
  1530, 1460, 1440, 1410, 1410, 1410, 1410
};
const G4double G4NucleonNuclearCrossSection::w_p_t[30] =   
{ 
  480,  900, 1500, 2350, 3020, 3420, 3650, 3775, 3875, 3830, 3750, 3700, 3630, 3550, 
  3550, 3290, 3070, 2890, 2840, 2730, 2725, 2720, 2770, 2805, 2828, 2865, 2770, 2640, 
  2570, 2490
};
const G4double G4NucleonNuclearCrossSection::w_p_in[30] =  
{ 
  325,  680,  990, 1500, 1850, 2150, 2250, 2300, 2350, 2330, 2280, 2230, 2200, 2120, 
  2130, 1900, 1780, 1670, 1635, 1600, 1602, 1605, 1610, 1615, 1630, 1660, 1620, 1550, 
  1530, 1460
};
const G4double G4NucleonNuclearCrossSection::e7[35] =      
{
  0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.12, 0.14, 0.16, 0.18,  0.2, 0.22, 
  0.25,  0.3, 0.35,  0.4, 0.45,  0.5, 0.55,  0.6,  0.7,  0.8,  0.9,    1,    2,    3,    
  5,   10,   20,   50,  100,  500, 1000
};
const G4double G4NucleonNuclearCrossSection::pb_m_t[35] =  
{
  5890, 5700, 5610, 5580, 5550, 5480, 5400, 5300, 5100, 4930, 4750, 4600, 4400, 4280, 
  4170, 3915, 3650, 3470, 3260, 3150, 3120, 3070, 3085, 3100, 3120, 3160, 3070, 2930, 
  2820, 2750, 2710, 2655, 2640, 2640, 2640
};
const G4double G4NucleonNuclearCrossSection::pb_m_in[35] = 
{
  1575, 2025, 2300, 2575, 2850, 3000, 3115, 3180, 3080, 2940, 2800, 2670, 2550, 2450, 
  2370, 2220, 2110, 2000, 1920, 1880, 1850, 1800, 1805, 1810, 1820, 1840, 1800, 1720, 
  1640, 1620, 1570, 1530, 1530, 1530, 1530
};
const G4double G4NucleonNuclearCrossSection::pb_p_t[30] =  
{ 
  515,  940, 1500, 2400, 3270, 3750, 4050, 4140, 4260, 4200, 4080, 3990, 3990, 3810, 
  3730, 3520, 3370, 3186, 3110, 3010, 2990, 2985, 3005, 3020, 3040, 3080, 3020, 2905, 
  2790, 2750
};
const G4double G4NucleonNuclearCrossSection::pb_p_in[30] = 
{ 
  348,  707, 1040, 1650, 2100, 2400, 2580, 2640, 2650, 2520, 2410, 2300, 2250, 2190, 
  2130, 2000, 1930, 1870, 1830, 1790, 1770, 1765, 1775, 1780, 1790, 1800, 1775, 1710, 
  1620, 1620
};
const G4double G4NucleonNuclearCrossSection::u_m_t[35] =   
{
  7080, 6830, 6650, 6530, 6400, 6280, 6100, 5840, 5660, 5520, 5330, 5160, 4990, 4810, 
  4630, 4323, 4130, 3870, 3700, 3550, 3490, 3465, 3467, 3475, 3495, 3515, 3440, 3360, 
  3150, 3040, 2985, 2955, 2940, 2940, 2940
};
const G4double G4NucleonNuclearCrossSection::u_m_in[35] =  
{
  1740, 2220, 2500, 2820, 3080, 3300, 3420, 3500, 3420, 3330, 3200, 3060, 2940, 2850, 
  2710, 2470, 2380, 2250, 2160, 2080, 2040, 2045, 2047, 2050, 2055, 2060, 2010, 1980, 
  1830, 1780, 1735, 1710, 1700, 1700, 1700
};
const G4double G4NucleonNuclearCrossSection::u_p_t[30] =   
{ 
  485,  960, 1580, 2700, 3550, 4050, 4320, 4420, 4620, 4660, 4580, 4470, 4350, 4295, 
  4187, 3938, 3755, 3573, 3450, 3342, 3310, 3295, 3310, 3330, 3375, 3405, 3350, 3338, 
  3135, 3040
};
const G4double G4NucleonNuclearCrossSection::u_p_in[30] =  
{ 
  334,  720, 1020, 1560, 2100, 2300, 2550, 2700, 2880, 2880, 2760, 2660, 2550, 2510, 
  2430, 2270, 2130, 2060, 2000, 1970, 1950, 1950, 1960, 1960, 1970, 1980, 1950, 1978, 
  1830, 1780
};

//////////////////////////////////////////////////////////////////////////////////
//
//

G4NucleonNuclearCrossSection::G4NucleonNuclearCrossSection()
{
   thePimData.push_back(new G4PiData(he_t,   he_in,  e1, 38));
   thePipData.push_back(new G4PiData(he_t,   he_in,  e1, 38));
   thePimData.push_back(new G4PiData(be_m_t, be_m_in, e1, 38));
   thePipData.push_back(new G4PiData(be_p_t, be_p_in, e1, 24));
   thePimData.push_back(new G4PiData(c_m_t,  c_m_in,  e2, 39));
   thePipData.push_back(new G4PiData(c_p_t,  c_p_in,  e2, 24));
   thePimData.push_back(new G4PiData(n_m_t,  n_m_in,  e2, 39));
   thePipData.push_back(new G4PiData(n_p_t,  n_p_in,  e2, 27));
   thePimData.push_back(new G4PiData(o_m_t,  o_m_in,  e3, 31));
   thePipData.push_back(new G4PiData(o_p_t,  o_p_in,  e3, 20));
   thePimData.push_back(new G4PiData(na_m_t, na_m_in, e3, 31));
   thePipData.push_back(new G4PiData(na_p_t, na_p_in, e3, 22));
   thePimData.push_back(new G4PiData(al_m_t, al_m_in, e3_1, 31));
   thePipData.push_back(new G4PiData(al_p_t, al_p_in, e3_1, 21));
   thePimData.push_back(new G4PiData(ca_m_t, ca_m_in, e3_1, 31));
   thePipData.push_back(new G4PiData(ca_p_t, ca_p_in, e3_1, 23));
   thePimData.push_back(new G4PiData(fe_m_t, fe_m_in, e4, 32));
   thePipData.push_back(new G4PiData(fe_p_t, fe_p_in, e4, 25));
   thePimData.push_back(new G4PiData(cu_m_t, cu_m_in, e4, 32));
   thePipData.push_back(new G4PiData(cu_p_t, cu_p_in, e4, 25));
   thePimData.push_back(new G4PiData(mo_m_t, mo_m_in, e5, 34));
   thePipData.push_back(new G4PiData(mo_p_t, mo_p_in, e5, 27));
   thePimData.push_back(new G4PiData(cd_m_t, cd_m_in, e5, 34));
   thePipData.push_back(new G4PiData(cd_p_t, cd_p_in, e5, 28));
   thePimData.push_back(new G4PiData(sn_m_t, sn_m_in, e6, 35));
   thePipData.push_back(new G4PiData(sn_p_t, sn_p_in, e6, 29));
   thePimData.push_back(new G4PiData(w_m_t,  w_m_in,  e6, 35));
   thePipData.push_back(new G4PiData(w_p_t,  w_p_in,  e6, 30));
   thePimData.push_back(new G4PiData(pb_m_t, pb_m_in, e7, 35));
   thePipData.push_back(new G4PiData(pb_p_t, pb_p_in, e7, 30));
   thePimData.push_back(new G4PiData(u_m_t,  u_m_in,  e7, 35));
   thePipData.push_back(new G4PiData(u_p_t,  u_p_in,  e7, 30));

   theZ.push_back(2); // He
   theZ.push_back(4); // Be
   theZ.push_back(6); // C
   theZ.push_back(7); // N
   theZ.push_back(8); // O
   theZ.push_back(11); // Na
   theZ.push_back(13); // Al
   theZ.push_back(20); // Ca
   theZ.push_back(26); // Fe
   theZ.push_back(29); // Cu
   theZ.push_back(42); // Mo
   theZ.push_back(48); // Cd
   theZ.push_back(50); // Sn
   theZ.push_back(74); // W
   theZ.push_back(82); // Pb
   theZ.push_back(92); // U

}

///////////////////////////////////////////////////////////////////////////////////
//
//

G4NucleonNuclearCrossSection::~G4NucleonNuclearCrossSection()
{
   std::for_each(thePimData.begin(), thePimData.end(), G4PiData::Delete());
   std::for_each(thePipData.begin(), thePipData.end(), G4PiData::Delete());
}

////////////////////////////////////////////////////////////////////////////
//
//

G4double G4NucleonNuclearCrossSection::
GetCrossSection( const G4DynamicParticle* aParticle, 
                 const G4Element* anElement,
                       G4double                )
{
   // precondition
   G4ping debug("debug_NucleonNuclearCrossSection");
	 G4double tuning = 1.;
   using namespace std;
   G4bool ok = false;
   if(aParticle->GetDefinition() == G4PionMinus::PionMinus()) ok=true;
   if(aParticle->GetDefinition() == G4PionPlus::PionPlus())   ok=true;
   if(!ok) 
   {
     throw G4HadronicException(__FILE__, __LINE__,
                          	"Call to G4NucleonNuclearCrossSection failed.");
   }
   G4double charge = aParticle->GetDefinition()->GetPDGCharge();
   G4double kineticEnergy = aParticle->GetKineticEnergy();

   // body
   G4double result = 0;
   G4int Z=G4lrint(anElement->GetZ());
   debug.push_back(Z);
   size_t it=0;
   while(it<theZ.size() && Z>theZ[it]) it++;
   debug.push_back(theZ[it]);
   debug.push_back(kineticEnergy);
	 if(Z==29) tuning=.96;
   if(Z > theZ[it]) 
   {
     throw G4HadronicException(__FILE__, __LINE__,
       "Called G4NucleonNuclearCrossSection outside parametrization");
   }
   G4int Z1, Z2;
   G4double x1, x2;
   if(charge<0)
   {
     if(theZ[it]==Z)
     {
       result = thePimData[it]->ReactionXSection(kineticEnergy);
       debug.push_back("D1 ");
       debug.push_back(result);
     }
     else
     {
       x1 = thePimData[it-1]->ReactionXSection(kineticEnergy);
       Z1 = theZ[it-1];
       x2 = thePimData[it]->ReactionXSection(kineticEnergy);
       Z2 = theZ[it];
       result = Interpolate(Z1, Z2, Z, x1, x2);
       debug.push_back("D2 ");
       debug.push_back(x1);
       debug.push_back(x2);
       debug.push_back(Z1);
       debug.push_back(Z2);
       debug.push_back(result);
     }
   }
   else
   {
     if(theZ[it]==Z)
     {
                              // at high energies, when no data for pi+, use pi- 
       std::vector<G4PiData *> * theData = &thePimData;
       if(thePipData[it]->AppliesTo(kineticEnergy))
       {
         theData = &thePipData;
       }
       result = theData->operator[](it)->ReactionXSection(kineticEnergy);
       debug.push_back("D3 ");
       debug.push_back(result);
     }
     else
     {
       std::vector<G4PiData *> * theLData = &thePimData;
       if(thePipData[it-1]->AppliesTo(kineticEnergy))
       {
         theLData = &thePipData;
       }
       std::vector<G4PiData *> * theHData = &thePimData;
       if(thePipData[it]->AppliesTo(kineticEnergy))
       {
         theHData = &thePipData;
       }
       x1 = theLData->operator[](it-1)->ReactionXSection(kineticEnergy);
       Z1 = theZ[it-1];
       x2 = theHData->operator[](it)->ReactionXSection(kineticEnergy);
       Z2 = theZ[it];
       result = Interpolate(Z1, Z2, Z, x1, x2);
       debug.push_back("D4 ");
       debug.push_back(x1);
       debug.push_back(x2);
       debug.push_back(Z1);
       debug.push_back(Z2);
       debug.push_back(result);
     }
   }

   debug.dump();
   return tuning*result;
}

/////////////////////////////////////////////////////////////////////////////
//
//

G4double G4NucleonNuclearCrossSection::
 Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2)
{ 
//   Nucleon numbers obtained from G4NistManager G4 8.0
 static const G4double A[92] = 
{
  1.0001, 4.0000, 6.9241, 9.0000, 10.801, 12.011, 14.004, 16.004, 19.000, 20.188,
  23.000, 24.320, 27.000, 28.109, 31.000, 32.094, 35.484, 39.985, 39.135, 40.116,
  45.000, 47.918, 50.998, 52.055, 55.000, 55.910, 59.000, 58.760, 63.617, 65.468,
  69.798, 72.691, 75.000, 79.042, 79.986, 83.887, 85.557, 87.710, 89.000, 91.318,
  93.000, 96.025, 98.000, 101.16, 103.00, 106.51, 107.96, 112.51, 114.91, 118.81,
  121.86, 127.70, 127.00, 131.39, 133.00, 137.42, 139.00, 140.21, 141.00, 144.32,
  145.00, 150.45, 152.04, 157.33, 159.00, 162.57, 165.00, 167.32, 169.00, 173.10,
  175.03, 178.54, 181.00, 183.89, 186.25, 190.27, 192.25, 195.11, 197.00, 200.63,
  204.41, 207.24, 209.00, 209.00, 210.00, 222.00, 223.00, 226.00, 227.00, 232.00,
  231.00, 237.98
};
			 
  static G4bool NeedInit=true;		     
  static G4double A75[92];

  if ( NeedInit )
  {
    for (G4int i=0; i<92; ++i)
    {
       A75[i]=std::pow(A[i],0.75);
    }
    NeedInit=false;
  }

//        for tabulated data, cross section scales with A^.75
   G4double r1 = x1 / A75[Z1-1] * A75[Z-1];
   G4double r2 = x2 / A75[Z2-1] * A75[Z-1];
   G4double result=0.5*(r1+r2);
//       G4cout << "x1/2, z1/2 z" <<x1<<" "<<x2<<" "<<Z1<<" "<<Z2<<" "<<Z<<G4endl;
//       G4cout << "res1/2 " << r1 <<" " << r2 <<" " << result<< G4endl;
   return result;
}
