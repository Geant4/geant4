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
// Based on G4NucleonNuclearCrossSection class
//
//

#include "G4ComponentBarNucleonNucleusXsc.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Pow.hh"

// Group 1: He, Be, C for 44 energies  

const G4double G4ComponentBarNucleonNucleusXsc::e1[44] =     
{
  0.014, 0.015, 0.017, 0.02, 0.022, 0.025, 0.027, 0.03, 0.035, 0.04,
  0.045, 0.05,  0.06,  0.07, 0.08,  0.09,  0.1,   0.12, 0.14,  0.15,
  0.16,  0.18,  0.20,  0.25, 0.30,  0.35,  0.4,   0.5,  0.6,   0.7,
  0.8,   0.9,   1.0,   1.5,  2.0,   3.0,   5.0,   7.0, 10.0,  20.0,
 50.0, 100.0, 500.0, 1000.0
};

const G4double G4ComponentBarNucleonNucleusXsc::he_m_t[44] =   
{ 
  1090, 1020, 915, 800, 710, 640, 600, 560, 500, 440, 390, 360, 295, 256, 220, 192, 
  168, 136, 120, 116, 114, 110, 107, 104, 106, 108, 110, 120, 126, 135, 140, 144, 146, 
  148, 152, 150, 146, 142, 138, 132, 129, 126, 127, 128  
};
const G4double G4ComponentBarNucleonNucleusXsc::he_m_in[44] =  
{ 
  0, 5, 10, 20, 35, 55, 70, 80, 90, 105, 115, 115, 100, 90, 86, 84, 84, 82, 80, 80, 80, 80, 
  79, 78, 80, 84, 88, 94, 100, 105, 108, 108, 108, 112, 114, 114, 112, 110, 108, 106, 104, 
  101, 102, 102
};
const G4double G4ComponentBarNucleonNucleusXsc::he_p_in[44] =  
{ 
  0, 2, 3, 13, 30, 50, 65, 77, 90, 105, 115, 115, 100, 90, 86, 84, 84, 82, 80, 80, 80, 80, 
  79, 78, 80, 84, 88, 94, 100, 105, 108, 108, 108, 112, 114, 114, 112, 110, 108, 106, 104, 
  101, 102, 102
};

const G4double G4ComponentBarNucleonNucleusXsc::be_m_t[44] = 
{
  1490, 1460, 1400, 1350, 1270, 1200, 1160, 1100, 1000, 910, 810, 740, 625, 575, 455, 406, 
  365, 310, 275, 262, 255, 240, 235, 225, 225, 230, 238, 252, 270, 282, 288, 290, 294, 303, 
  303, 300, 292, 284, 277, 267, 263, 264, 268, 268 
};
const G4double G4ComponentBarNucleonNucleusXsc::be_m_in[44] =
{ 
  650, 640, 617, 595, 555, 520, 495, 470, 430, 385, 350, 320, 270, 250, 210, 190, 185, 178, 
  175, 175, 175, 175, 175, 170, 170, 172, 176, 184, 194, 200, 209, 213, 214, 216, 216, 212, 
  210, 210, 210, 210, 210, 210, 210, 210
};
const G4double G4ComponentBarNucleonNucleusXsc::be_p_in[44] =
{ 
  490, 540, 580, 545, 525, 495, 470, 450, 420, 370, 340, 310, 262, 242, 205, 185, 180, 175, 
  172, 175, 175, 175, 175, 170, 170, 172, 176, 184, 194, 200, 209, 213, 214, 216, 216, 212, 
  210, 210, 210, 210, 210, 210, 210, 210
};

const G4double G4ComponentBarNucleonNucleusXsc::c_m_t[44] = 
{
  1240, 1370, 1450, 1455, 1445, 1385, 1345, 1290, 1210, 1110, 1020, 940, 800, 700, 604, 530, 
  475, 396, 350, 336, 320, 303, 294, 280, 280, 286, 296, 314, 330, 344, 356, 360, 364, 384, 
  388, 384, 364, 352, 344, 330, 324, 324, 332, 332
};
const G4double G4ComponentBarNucleonNucleusXsc::c_m_in[44] =
{
  590, 570, 542, 510, 500, 460, 445, 430, 395, 380, 350, 330, 295, 270, 255, 240, 228, 222, 
  216, 216, 210, 210, 210, 208, 210, 214, 216, 228, 240, 248, 254, 257, 260, 262, 260, 256, 
  252, 252, 250, 250, 248, 248, 248, 248
};
const G4double G4ComponentBarNucleonNucleusXsc::c_p_in[44] =
{ 
  310, 330, 400, 440, 450, 435, 430, 420, 385, 370, 340, 320, 288, 263, 249, 234, 222, 216, 
  210, 211, 205, 208, 210, 208, 210, 214, 216, 228, 240, 248, 254, 257, 260, 262, 260, 256, 
  252, 252, 250, 250, 248, 248, 248, 248
};

// Group 2: N, O, Na for 44 energies (e1=e2)

const G4double G4ComponentBarNucleonNucleusXsc::e2[44] =    
{
 0.014, 0.015, 0.017, .02, 0.022, 0.025, 0.027, 0.03, 0.035, .04, 0.045, 0.05, .06, 0.07, 
  .08, 0.09,  .1, .12, .14, .15, .16, .18, .20, .25, .30, .35, .4 , 0.5, 0.6, 0.7,  0.8,  
  0.9,   1, 1.5,   2,   3,   5,  7, 10,   
  20,   50,  100,  500, 1000  
};

const G4double G4ComponentBarNucleonNucleusXsc::n_m_t[44] = 
{
  1420,1480, 1537, 1550, 1525, 1500, 1480, 1425, 1340, 1260, 1175, 1090, 930, 805, 690, 612, 
  552, 462, 402, 384, 372, 350, 345, 326, 324, 328, 336, 356, 372, 388, 400, 408, 415, 430, 
  435, 432, 415, 402, 390, 375, 367, 370, 382, 385
};
const G4double G4ComponentBarNucleonNucleusXsc::n_m_in[44] =
{
  680, 665, 625, 580, 562, 525, 510, 485, 450, 435, 410, 387, 340, 310, 290, 280, 276, 274, 
  260, 258, 254, 247, 245, 240, 240, 244, 250, 260, 268, 275, 280, 285, 290, 295, 300, 294, 
  292, 290, 285, 285, 282, 282, 282, 282
};
const G4double G4ComponentBarNucleonNucleusXsc::n_p_in[44] =
{ 
  420, 440, 470, 490, 497, 500, 480, 462, 440, 425, 400, 377, 333, 303, 284, 274, 270, 268, 
  254, 252, 247, 245, 245, 240, 240, 244, 250, 260, 268, 275, 280, 285, 290, 295, 300, 294, 
  292, 290, 285, 285, 282, 282, 282, 282
};

const G4double G4ComponentBarNucleonNucleusXsc::o_m_t[44] =  
{
  1520, 1570, 1630, 1660, 1647, 1623, 1595, 1555, 1475, 1395, 1290, 1207, 1035, 925, 816, 
  720, 645, 540, 462, 438, 415, 392, 378, 362, 361, 381, 390, 403, 417, 440, 460, 470, 
  479, 498, 504, 498, 477, 457, 443, 427, 420, 425, 429, 430
};
const G4double G4ComponentBarNucleonNucleusXsc::o_m_in[44] = 
{
  750, 740, 700, 650, 620, 575, 555, 530, 505, 462, 435, 420, 375, 345, 320, 310, 300, 293, 
  288, 282, 282, 280, 276, 270, 271, 275, 280, 290, 295, 304, 310, 315, 318, 332, 335, 330, 
  323, 320, 317, 315, 315, 315, 315, 315
};
const G4double G4ComponentBarNucleonNucleusXsc::o_p_in[44] = 
{
  460, 485, 510, 535, 537, 532, 520, 500, 460, 432, 405, 390, 350, 320, 310, 304, 293, 287, 
  283, 279, 279, 278, 276, 270, 271, 275, 280, 290, 295, 304, 310, 315, 318, 332, 335, 330, 
  323, 320, 317, 315, 315, 315, 315, 315
};

const G4double G4ComponentBarNucleonNucleusXsc::na_m_t[44] = 
{
  1570, 1620, 1695, 1730, 1750, 1760, 1755, 1740, 1710, 1643, 1560, 1480, 1343, 1220, 1073, 
  953, 860, 720, 618, 582, 546, 522, 504, 484, 492, 500, 512, 538, 560, 586, 608, 622, 632, 
  660, 668, 664, 640, 616, 596, 568, 568, 568, 568, 568 
};
const G4double G4ComponentBarNucleonNucleusXsc::na_m_in[44] =
{
  960, 930, 890, 822, 790, 750, 725, 686, 620, 600, 575, 540, 497, 450, 414, 390, 380, 372, 
  354, 360, 355, 354, 350, 350, 350, 356, 364, 384, 392, 400, 408, 410, 420, 408, 412, 420, 
  411, 409, 407, 403, 400, 400, 400, 400
};
const G4double G4ComponentBarNucleonNucleusXsc::na_p_in[44] =
{
  600, 617, 660, 675, 680, 680, 670, 650, 575, 550, 525, 490, 450, 420, 385, 367, 360, 350, 
  350, 350, 345, 347, 350, 350, 350, 356, 364, 384, 392, 400, 408, 410, 420, 408, 412, 420, 
  411, 409, 407, 403, 400, 400, 400, 400
};

// Al, Si, Ca for 45 energies

const G4double G4ComponentBarNucleonNucleusXsc::e3[45] =     
{
  0.014, 0.015, 0.016, 0.017, .02, 0.022, 0.025, 0.027, 0.03, 0.035, .04, 0.045, 0.05, .06, 0.07, 
    .08, 0.09,   .1,    .12,  .14,  .15,   .16,   .18,   .20,  .25,  .30,  .35,  0.4, 0.5,  0.6, 
    0.7, 0.8,   0.9,    1,   1.5,   2,     3,     5,     7,   10,   20,   50,  100,   500, 1000  
};

const G4double G4ComponentBarNucleonNucleusXsc::al_m_t[45] = 
{ 
  1735, 1750, 1760, 1795, 1830, 1855, 1885, 1895, 1900, 1870, 1835, 1785, 1710, 1522, 1350, 
  1212, 1080,  972,  816,  720,  678,  642,  600,  567,  558,  560,  578,  592,  616,  644,  
   672,  688,  708,  720,  736,  754,  736,  706,  680,  672,  646,  632,  632,  632,  632
};
const G4double G4ComponentBarNucleonNucleusXsc::al_m_in[45] = 
{
  1000,  990,  975,  950,  905,  875,  825,  800,  762,  690,  652,  610,  570,  495,  480, 
   456,  444,  432,  420,  420,  420,  420,  410,  410,  400,  402,  404,  408,  424,  438, 
   448,  450,  454,  456,  472,  480,  466,  456,  452,  448,  444,  440,  440,  440,  440
};
const G4double G4ComponentBarNucleonNucleusXsc::al_p_in[45] = 
{
   650,  682,  690,  715,  750,  762,  750,  740,  720,  655,  617,  575,  540,  470,  455, 
   //   532,  420,  408,  400,  403,  403,  408,  406,  404,  400,  402,  404,  408,  424,  438, 
   432,  420,  408,  400,  403,  403,  408,  406,  404,  400,  402,  404,  408,  424,  438, 
   448,  450,  454,  456,  472,  480,  466,  456,  452,  448,  444,  440,  440,  440,  440 
};

const G4double G4ComponentBarNucleonNucleusXsc::si_m_t[45] = 
{ 
  1810, 1833, 1850, 1872, 1920, 1950, 1995, 2020, 2035, 2000, 1930, 1850, 1760, 1570, 1400, 
  1255, 1110, 1008,  846,  742,  696,  671,  623,  588,  584,  584,  602,  618,  645,  679, 
   708,  727,  746,  757,  769,  782,  771,  734,  710,  698,  672,  654,  650,  650,  650 
};
const G4double G4ComponentBarNucleonNucleusXsc::si_m_in[45] = 
{
  1060, 1035, 1015,  990,  935,  900,  860,  830,  790,  725,  665,  630,  600,  520,  504, 
   486,  470,  456,  444,  432,  432,  432,  418,  418,  415,  412,  416,  422,  440,  460, 
   472,  476,  479,  480,  492,  496,  488,  472,  472,  464,  460,  452,  448,  448,  448 
};
const G4double G4ComponentBarNucleonNucleusXsc::si_p_in[45] = 
{
   670,  700,  725,  750,  780,  780,  770,  757,  735,  690,  635,  585,  570,  490,  475, 
   460,  446,  431,  423,  425,  425,  425,  425,  422,  422,  412,  416,  422,  440,  460, 
   472,  476,  479,  480,  492,  496,  488,  472,  472,  464,  460,  452,  448,  448,  448
};

const G4double G4ComponentBarNucleonNucleusXsc::ca_m_t[45] = 
{ 
  2180, 2130, 2095, 2075, 2115, 2150, 2220, 2250, 2300, 2365, 2360, 2280, 2180, 2000, 
  1805, 1650, 1500, 1340, 1140, 990, 940, 890, 825, 790, 770, 773, 787, 800, 830, 870, 
  905, 930, 950, 965, 990, 1002, 990, 965, 945, 925, 892, 860, 860, 860, 860 
};
const G4double G4ComponentBarNucleonNucleusXsc::ca_m_in[45] = 
{
  1240, 1225, 1200, 1180, 1125, 1090, 1045, 1020, 980, 925, 880, 825, 770, 680, 640, 
  620, 615, 600, 580, 565, 560, 560, 560, 550, 535, 530, 540, 550, 570, 595, 610, 615, 
  620, 622, 629, 630, 620, 612, 607, 592, 587, 580, 580, 580, 580 
};
const G4double G4ComponentBarNucleonNucleusXsc::ca_p_in[45] = 
{
  770, 800, 823, 850, 900, 925, 935, 920, 895, 835, 800, 750, 715, 640, 605, 590, 588, 
  573, 555, 543, 540, 540, 540, 535, 530, 530, 540, 550, 570, 595, 610, 615, 
  620, 622, 629, 630, 620, 612, 607, 592, 587, 580, 580, 580, 580 
};

// Fe, Cu, Mo for 47 energies

const G4double G4ComponentBarNucleonNucleusXsc::e4[47] =     
{
  0.014, 0.015, 0.017, .02, 0.022, 0.025, 0.027, 0.03, 0.033, 0.035, 0.037, .04, 0.045, 
  0.05, 0.055, .06, 0.07, .08, 0.09,  .1, .12, .14, .15, .16, .18, .20, .25, .30, .35, 
  .4 , 0.5, 0.6, 0.7,  0.8,  0.9,   1, 1.5,   2,   3,   5,  7, 10,   
  20,   50,  100,  500, 1000
};

const G4double G4ComponentBarNucleonNucleusXsc::fe_m_t[47] = 
{
  2580, 2490, 2370, 2282, 2275, 2285, 2320, 2370, 2432, 2445, 2460, 2485, 2530, 2540, 
  2517, 2480, 2290, 2110, 1940, 1790, 1510, 1290, 1220, 1150, 1070, 1030, 1013, 1020, 
  1030, 1043, 1075, 1110, 1133, 1163, 1185, 1225, 1252, 1260, 1260, 1233, 1207, 1185, 
  1140, 1110, 1110, 1110, 1110
};
const G4double G4ComponentBarNucleonNucleusXsc::fe_m_in[47] = 
{
  1440, 1433, 1390, 1325, 1280, 1260, 1215, 1180, 1140, 1110, 1080, 1040, 990, 955, 920, 
  885, 835, 800, 780, 765, 750, 725, 720, 720, 710, 700, 700, 700, 712, 705, 735, 750, 
  765, 775, 780, 795, 810, 813, 810, 784, 757, 743, 735, 720, 720, 720, 720  
};
const G4double G4ComponentBarNucleonNucleusXsc::fe_p_in[47] = 
{
  900, 960, 1070, 1090, 1115, 1120, 1115, 1080, 1045, 1025, 1000, 960, 900, 885, 865, 790, 
  765, 740, 720, 700, 697, 697, 697, 697, 695, 690, 688, 690, 712, 705, 735, 750, 
  765, 775, 780, 795, 810, 813, 810, 784, 757, 743, 735, 720, 720, 720, 720 
};

const G4double G4ComponentBarNucleonNucleusXsc::cu_m_t[47] = 
{
  2920, 2800, 2615, 2480, 2455, 2430, 2440, 2460, 2500, 2530, 2560, 2615, 2690, 2720, 
  2700, 2645, 2500, 2320, 2140, 1970, 1670, 1460, 1380, 1285, 1200, 1160, 1140, 1147, 
  1163, 1170, 1200, 1237, 1265, 1285, 1305, 1328, 1375, 1390, 1395, 1370, 1335, 1315, 
  1270, 1230, 1230, 1230, 1230 
};
const G4double G4ComponentBarNucleonNucleusXsc::cu_m_in[47] = 
{
  1540, 1535, 1500, 1445, 1407, 1380, 1330, 1300, 1285, 1270, 1240, 1190, 1090, 1010, 
  940, 920, 860, 835, 820, 810, 800, 780, 775, 770, 760, 760, 758, 765, 765, 770, 795, 
  810, 825, 830, 840, 848, 870, 870, 868, 840, 825, 810, 803, 795, 795, 795, 795 
};
const G4double G4ComponentBarNucleonNucleusXsc::cu_p_in[47] = 
{
  935, 1000, 1060, 1190, 1220, 1250, 1240, 1210, 1150, 1130, 1115, 1050, 985, 950, 890, 
  870, 820, 800, 785, 780, 770, 750, 745, 740, 735, 735, 745, 760, 762, 770, 795, 
  810, 825, 830, 840, 848, 870, 870, 868, 840, 825, 810, 803, 795, 795, 795, 795 
};

const G4double G4ComponentBarNucleonNucleusXsc::mo_m_t[47] = 
{
  4150, 4040, 3800, 3490, 3300, 3060, 2960, 2845, 2785, 2820, 2850, 2980, 3170, 3230, 
  3270, 3280, 3225, 3075, 2895, 2710, 2355, 2060, 1925, 1800, 1630, 1560, 1540, 1550, 
  1570, 1590, 1650, 1685, 1715, 1740, 1760, 1780, 1850, 1880, 1858, 1815, 1790, 1782, 
  1720, 1690, 1690, 1690, 1690 
};
const G4double G4ComponentBarNucleonNucleusXsc::mo_m_in[47] = 
{
  1790, 1775, 1740, 1680, 1640, 1580, 1550, 1510, 1460, 1440, 1418, 1380, 1330, 1280, 
  1240, 1200, 1155, 1140, 1110, 1110, 1080, 1065, 1050, 1050, 1025, 1020, 1015, 1020, 
  1022, 1026, 1060, 1085, 1100, 1110, 1120, 1127, 1150, 1160, 1140, 1100, 1085, 1080, 
  1070, 1070, 1070, 1070, 1070  
};
const G4double G4ComponentBarNucleonNucleusXsc::mo_p_in[47] = 
{
  1025, 1080, 1190, 1380, 1440, 1495, 1475, 1420, 1350, 1310, 1300, 1290, 1250, 1200, 
  1170, 1130, 1095, 1060, 1040, 1022, 1020, 1016, 1016, 1016, 1016, 1012, 1005, 1005, 
  1005, 1010, 1060, 1085, 1100, 1110, 1120, 1127, 1150, 1160, 1140, 1100, 1085, 1080, 
  1070, 1070, 1070, 1070, 1070
};

// Cd, Sn, W for 48 energies

const G4double G4ComponentBarNucleonNucleusXsc::e5[48] =     
{
   0.014, 0.015, 0.017, 0.018, .02, 0.022, 0.025, 0.027, 0.03, 0.033, 0.035, .04, 
   0.045, 0.05, 0.055, .06, .065, 0.07, .08, 0.09,  .1, .12, .14, .15, .16, .18, 
   .20, .25, .30, .35, .4 , 0.5, 0.6, 0.7,  0.8, 0.9, 1, 1.5, 2, 3, 5, 7, 10,  20, 
   50, 100, 500, 1000
};

const G4double G4ComponentBarNucleonNucleusXsc::cd_m_t[48] = 
{
  4420, 4280, 4170, 4070, 3860, 3680, 3420, 3280, 3125, 3060, 3080, 3190, 3350, 3445, 
  3510, 3540, 3560, 3550, 3460, 3300, 3030, 2640, 2340, 2190, 2070, 1950, 1770, 1732, 
  1740, 1760, 1780, 1832, 1885, 1925, 1945, 1960, 1980, 2070, 2080, 2065, 2040, 2022, 
  1980, 1940, 1870, 1870, 1870, 1870 
};
const G4double G4ComponentBarNucleonNucleusXsc::cd_m_in[48]= 
{
  1920, 1910, 1880, 1860, 1840, 1800, 1760, 1720, 1675, 1630, 1600, 1520, 1465, 1420, 
  1390, 1340, 1310, 1280, 1275, 1235, 1225, 1200, 1170, 1170, 1170, 1165, 1145, 1140, 
  1140, 1135, 1160, 1180, 1220, 1240, 1250, 1260, 1265, 1270, 1275, 1250, 1222, 1222, 
  1220, 1215, 1190, 1190, 1190, 1190 
};
const G4double G4ComponentBarNucleonNucleusXsc::cd_p_in[48] = 
{
  1020, 1100, 1225, 1290, 1440, 1520, 1575, 1560, 1518, 1460, 1420, 1400, 1365, 1340, 
  1300, 1280, 1260, 1200, 1190, 1160, 1125, 1125, 1125, 1125, 1125, 1125, 1120, 1120, 
  1120, 1118, 1146, 1180, 1220, 1240, 1250, 1260, 1265, 1270, 1275, 1250, 1222, 1222, 
  1220, 1215, 1190, 1190, 1190, 1190 
};

const G4double G4ComponentBarNucleonNucleusXsc::sn_m_t[48] =  
{
  4420, 4400, 4260, 4150, 3980, 3770, 3530, 3370, 3245, 3180, 3170, 3260, 3400, 3500, 
  3560, 3610, 3650, 3680, 3580, 3390, 3190, 2760, 2430, 2295, 2175, 1990, 1880, 1810, 
  1820, 1840, 1865, 1940, 1985, 2020, 2040, 2060, 2080, 2160, 2185, 2180, 2110, 2105, 
  2080, 2050, 1980, 1980, 1980, 1980 
};
const G4double G4ComponentBarNucleonNucleusXsc::sn_m_in[48] = 
{
  1945, 1940, 1905, 1890, 1860, 1830, 1780, 1755, 1717, 1680, 1645, 1570, 1500, 1455, 
  1410, 1370, 1340, 1320, 1290, 1285, 1260, 1240, 1235, 1212, 1200, 1200, 1200, 1190, 
  1190, 1200, 1210, 1240, 1270, 1285, 1300, 1300, 1310, 1320, 1320, 1290, 1240, 1240, 
  1240, 1240, 1240, 1240, 1240, 1240 
};
const G4double G4ComponentBarNucleonNucleusXsc::sn_p_in[48] = 
{ 
  1020, 1080, 1270, 1335, 1465, 1505, 1610, 1610, 1550, 1535, 1500, 1440, 1407, 1370, 
  1340, 1300, 1285, 1260, 1230, 1215, 1200, 1180, 1170, 1170, 1165, 1165, 1170, 1165, 
  1165, 1183, 1195, 1240, 1270, 1285, 1300, 1300, 1310, 1320, 1320, 1290, 1240, 1240, 
  1240, 1240, 1240, 1240, 1240, 1240
};

const G4double G4ComponentBarNucleonNucleusXsc::w_m_t[48] =   
{
  5320, 5430, 5480, 5450, 5330, 5190, 4960, 4790, 4550, 4340, 4200, 4070, 4000, 4030, 
  4125, 4220, 4270, 4390, 4440, 4360, 4200, 3800, 3380, 3200, 3040, 2790, 2660, 2575, 
  2575, 2600, 2640, 2690, 2755, 2790, 2812, 2837, 2850, 2950, 3000, 2970, 2940, 2910, 
  2880, 2820, 2730, 2730, 2730, 2730 
};
const G4double G4ComponentBarNucleonNucleusXsc::w_m_in[48] =  
{
  2440, 2400, 2370, 2350, 2310, 2270, 2220, 2195, 2150, 2100, 2070, 2010, 1945, 1900, 
  1850, 1820, 1780, 1760, 1730, 1720, 1680, 1680, 1660, 1660, 1650, 1650, 1640, 1640, 
  1612, 1615, 1625, 1640, 1700, 1720, 1730, 1740, 1750, 1780, 1780, 1750, 1740, 1735, 
  1710, 1695, 1680, 1680, 1680, 1680  
};
const G4double G4ComponentBarNucleonNucleusXsc::w_p_in[48] =  
{ 
  950,  1020, 1240, 1400, 1560, 1670, 1760, 1830, 1850, 1855, 1870, 1840, 1800, 1770, 
  1740, 1715, 1680, 1670, 1650, 1620, 1610, 1600, 1600, 1600, 1600, 1600, 1600, 1595, 
  1585, 1595, 1615, 1640, 1700, 1720, 1730, 1740, 1750, 1780, 1780, 1750, 1740, 1735, 
  1710, 1695, 1680, 1680, 1680, 1680
};

// Pb, U for 46 energies

const G4double G4ComponentBarNucleonNucleusXsc::e6[46] =      
{
  0.014, 0.015, 0.017, 0.019, 0.02, 0.022, 0.025, 0.027, 0.03, 0.035,
  0.04,  0.045, 0.05,  0.055, 0.06, 0.07,  0.08,  0.09,  0.1,  0.12,
  0.14,  0.15,  0.16,  0.18,  0.20, 0.25,  0.30,  0.35,  0.4 , 0.5,
  0.6,   0.7,   0.8,   0.9,   1.0,  1.5,   2.0,   3.0,   5.0,  7.0,
 10.0,  20.0,  50.0, 100.0, 500.0, 1000.0
};

const G4double G4ComponentBarNucleonNucleusXsc::pb_m_t[46] =  
{
  5300, 5440, 5720, 5880, 5765, 5745, 5480, 5280, 4970, 4550, 4390, 4300, 4265, 4325, 
  4450, 4540, 4740, 4710, 4600, 4100, 3660, 3480, 3300, 3000, 2890, 2865, 2855, 2850, 
  2865, 2920, 2955, 3000, 3030, 3060, 3105, 3240, 3290, 3270, 3240, 3180, 3090, 3060, 
  2970, 2970, 2970, 2970  

};
const G4double G4ComponentBarNucleonNucleusXsc::pb_m_in[46] = 
{
  2580, 2550, 2505, 2462, 2460, 2435, 2380, 2355, 2280, 2180, 2170, 2130, 2080, 2035, 
  1980, 1940, 1900, 1870, 1840, 1800, 1800, 1800, 1780, 1760, 1760, 1740, 1730, 1725, 
  1740, 1785, 1815, 1835, 1860, 1890, 1895, 1920, 1920, 1890, 1850, 1835, 1830, 1830, 
  1830, 1830, 1830, 1830 
};
const G4double G4ComponentBarNucleonNucleusXsc::pb_p_in[46] = 
{ 
  900,  1060, 1200, 1420, 1515, 1620, 1750, 1800, 1915, 2030, 1960, 1940, 1910, 1860, 
  1840, 1780, 1770, 1760, 1740, 1720, 1725, 1740, 1740, 1730, 1720, 1700, 1710, 1720, 
  1730, 1740, 1815, 1835, 1860, 1890, 1895, 1920, 1920, 1890, 1850, 1835, 1830, 1830, 
  1830, 1830, 1830, 1830
};

const G4double G4ComponentBarNucleonNucleusXsc::u_m_t[46] =   
{
  5800, 5940, 6160, 6345, 6360, 6350, 6170, 6020, 5760, 5350, 4990, 4800, 4710, 4690, 
  4760, 5040, 5190, 5200, 5080, 4600, 4120, 3920, 3720, 3420, 3240, 3150, 3160, 3180, 
  3210, 3240, 3280, 3350, 3390, 3435, 3480, 3560, 3585, 3580, 3540, 3500, 3470, 3410, 
  3335, 3335, 3335, 3335   
};
const G4double G4ComponentBarNucleonNucleusXsc::u_m_in[46] =  
{
  2820, 2770, 2700, 2660, 2645, 2620, 2580, 2550, 2515, 2450, 2390, 2320, 2260, 2225, 
  2200, 2140, 2080, 2060, 2040, 2000, 1980, 1965, 1960, 1930, 1920, 1890, 1905, 1920, 
  1945, 1970, 1985, 2010, 2040, 2070, 2080, 2090, 2095, 2080, 2063, 2060, 2050, 2040, 
  2005, 2005, 2005, 2005 
};
const G4double G4ComponentBarNucleonNucleusXsc::u_p_in[46] =  
{ 
  800,  900,  1100, 1300, 1410, 1510, 1680, 1800, 2000, 2200, 2080, 2060, 2035, 2100, 
  2030, 2030, 2000, 1960, 1960, 1960, 1940, 1925, 1920, 1905, 1890, 1860, 1880, 1910, 
  1930, 1945, 1985, 2010, 2040, 2070, 2080, 2090, 2095, 2080, 2063, 2060, 2050, 2040, 
  2005, 2005, 2005, 2005
};

using namespace std;

///////////////////////////////////////////////////////////////////////////////

G4ComponentBarNucleonNucleusXsc::G4ComponentBarNucleonNucleusXsc()
 : G4VComponentCrossSection("G4ComponentBarNucleonNucleusXsc"),
   fTotalXsc(0.0), fInelasticXsc(0.0), fElasticXsc(0.0)
{
  theNeutron = G4Neutron::Neutron();
  theProton  = G4Proton::Proton();
  
  // He, Be, C
   
   thePimData.push_back(new G4PiData(he_m_t, he_m_in, e1, 44));
   thePipData.push_back(new G4PiData(he_m_t, he_p_in, e1, 44));

   thePimData.push_back(new G4PiData(be_m_t, be_m_in, e1, 44));
   thePipData.push_back(new G4PiData(be_m_t, be_p_in, e1, 44));

   thePimData.push_back(new G4PiData(c_m_t,  c_m_in,  e1, 44));
   thePipData.push_back(new G4PiData(c_m_t,  c_p_in,  e1, 44));

   // N, O, Na

   thePimData.push_back(new G4PiData(n_m_t,  n_m_in,  e2, 44));
   thePipData.push_back(new G4PiData(n_m_t,  n_p_in,  e2, 44));

   thePimData.push_back(new G4PiData(o_m_t,  o_m_in,  e2, 44));
   thePipData.push_back(new G4PiData(o_m_t,  o_p_in,  e2, 44));

   thePimData.push_back(new G4PiData(na_m_t, na_m_in, e2, 44));
   thePipData.push_back(new G4PiData(na_m_t, na_p_in, e2, 44));

   // Al, Si, Ca

   thePimData.push_back(new G4PiData(al_m_t, al_m_in, e3, 45));
   thePipData.push_back(new G4PiData(al_m_t, al_p_in, e3, 45));

   thePimData.push_back(new G4PiData(si_m_t, si_m_in, e3, 45));
   thePipData.push_back(new G4PiData(si_m_t, si_p_in, e3, 45));

   thePimData.push_back(new G4PiData(ca_m_t, ca_m_in, e3, 45));
   thePipData.push_back(new G4PiData(ca_m_t, ca_p_in, e3, 45));

   // Fe, Cu, Mo

   thePimData.push_back(new G4PiData(fe_m_t, fe_m_in, e4, 47));
   thePipData.push_back(new G4PiData(fe_m_t, fe_p_in, e4, 47));

   thePimData.push_back(new G4PiData(cu_m_t, cu_m_in, e4, 47));
   thePipData.push_back(new G4PiData(cu_m_t, cu_p_in, e4, 47));

   thePimData.push_back(new G4PiData(mo_m_t, mo_m_in, e4, 47));
   thePipData.push_back(new G4PiData(mo_m_t, mo_p_in, e4, 47));

   // Cd, Sn, W

   thePimData.push_back(new G4PiData(cd_m_t, cd_m_in, e5, 48));
   thePipData.push_back(new G4PiData(cd_m_t, cd_p_in, e5, 48));

   thePimData.push_back(new G4PiData(sn_m_t, sn_m_in, e5, 48));
   thePipData.push_back(new G4PiData(sn_m_t, sn_p_in, e5, 48));

   thePimData.push_back(new G4PiData(w_m_t,  w_m_in,  e5, 48));
   thePipData.push_back(new G4PiData(w_m_t,  w_p_in,  e5, 48));

   // Pb, U

   thePimData.push_back(new G4PiData(pb_m_t, pb_m_in, e6, 46));
   thePipData.push_back(new G4PiData(pb_m_t, pb_p_in, e6, 46));

   thePimData.push_back(new G4PiData(u_m_t,  u_m_in,  e6, 46));
   thePipData.push_back(new G4PiData(u_m_t,  u_p_in,  e6, 46));

   theZ.push_back(2); // He
   theZ.push_back(4); // Be
   theZ.push_back(6); // C
   theZ.push_back(7); // N
   theZ.push_back(8); // O
   theZ.push_back(11); // Na
   theZ.push_back(13); // Al
   theZ.push_back(14); // Si
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

///////////////////////////////////////////////////////////////////////////////
//

G4ComponentBarNucleonNucleusXsc::~G4ComponentBarNucleonNucleusXsc()
{
   std::for_each(thePimData.begin(), thePimData.end(), G4PiData::Delete());
   std::for_each(thePipData.begin(), thePipData.end(), G4PiData::Delete());
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fTotalXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int )
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double )
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fInelasticXsc;
}

//////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double )
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fElasticXsc;
}

///////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int )
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
						 kinEnergy);
  fInelasticXsc = GetElementCrossSection(aDP, Z);
  delete aDP;

  return fElasticXsc;
}



////////////////////////////////////////////////////////////////////////////
//

G4bool
G4ComponentBarNucleonNucleusXsc::IsElementApplicable(const G4DynamicParticle* aParticle, 
						     G4int Z) // , const G4Material*)
{
  G4bool result = false;
  if(aParticle->GetDefinition() == theNeutron ) result = true;
  if(aParticle->GetDefinition() == theProton)   result = true;
  if(Z < 2)                                     result = false;
  if(aParticle->GetKineticEnergy() > 999.9*GeV) result = false;
  return result;
}

////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentBarNucleonNucleusXsc::GetElementCrossSection(const G4DynamicParticle* aParticle, 
							G4int Z) // , const G4Material*)
{
   G4double kineticEnergy = aParticle->GetKineticEnergy();
  
   G4double result = 0;
   // G4cout<<"Z = "<<Z<<G4endl;

   size_t it = 0;
   size_t itmax = theZ.size() - 1;
   for(; it <= itmax; ++it) { if(Z <= theZ[it]) { break; } }
   if( it > itmax ) { it = itmax; }
   G4int Z1, Z2;
   G4double x1, x2, xt1, xt2;

   std::vector<G4PiData *> * theData = &thePimData;
   if(aParticle->GetDefinition() == theProton) { theData = &thePipData; }

   if( theZ[it] == Z )
     {
       result = (*theData)[it]->ReactionXSection(kineticEnergy);
       fTotalXsc = (*theData)[it]->TotalXSection(kineticEnergy);
     }
   else
     {
       if(0 == it) { it = 1; }
       x1  = (*theData)[it-1]->ReactionXSection(kineticEnergy);
       xt1 = (*theData)[it-1]->TotalXSection(kineticEnergy);
       Z1  = theZ[it-1];
       x2  = (*theData)[it]->ReactionXSection(kineticEnergy);
       xt2 = (*theData)[it]->TotalXSection(kineticEnergy);
       Z2  = theZ[it];

       result = Interpolate(Z1, Z2, Z, x1, x2);
       fTotalXsc = Interpolate(Z1, Z2, Z, xt1, xt2);
     }

   fElasticXsc = fTotalXsc - result;
   if( fElasticXsc < 0.) { fElasticXsc = 0.; }

   return result;
}

/////////////////////////////////////////////////////////////////////////////
//

G4double G4ComponentBarNucleonNucleusXsc::
Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2)
{ 
//   Nucleon numbers obtained from G4NistManager G4 8.0

  static const G4double alpha = 2./3.;

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
  static G4ThreadLocal G4bool NeedInit = true;
		     
  static G4ThreadLocal G4double A75[92];

  if ( NeedInit )
  {
    for (G4int i=0; i<92; ++i)
    {
      A75[i] = G4Pow::GetInstance()->powA(A[i], alpha); // interpolate by square ~ A^(2/3)
    }
    NeedInit=false;
  }

  // for tabulated data, cross section scales with A^(2/3)
  G4double r1 = x1 / A75[Z1-1] * A75[Z-1];
  G4double r2 = x2 / A75[Z2-1] * A75[Z-1];
  G4double result = 0.5*(r1+r2);

  // More precise average
  if(Z1 != Z2) {
    G4double alp1 = (A[Z-1] - A[Z1-1]);
    G4double alp2 = (A[Z2-1] - A[Z-1]);
    result = (r1*alp2 + r2*alp1)/(alp1 + alp2);
  }
  //       G4cout << "x1/2, z1/2 z" <<x1<<" "<<x2<<" "<<Z1<<" "<<Z2<<" "<<Z<<G4endl;
  //       G4cout << "res1/2 " << r1 <<" " << r2 <<" " << result<< G4endl;
  return result;
}

void
G4ComponentBarNucleonNucleusXsc::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ComponentBarNucleonNucleusXsc is a variant of the Barashenkov\n"
          << "cross section parameterization to be used of protons and\n"
          << "nucleons on targets heavier than hydrogen.  It is intended for\n"
          << "use as a cross section component and is currently used by\n"
          << "G4BGGNucleonInelasticXS.  It is valid for incident energies up\n"
          << "to 1 TeV.\n"; 
}

