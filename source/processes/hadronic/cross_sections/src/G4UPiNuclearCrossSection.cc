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
//
// Calculation of the total, elastic and inelastic cross-sections
// based on Barashenkov parametrisations of pion data
//
// 16.08.06 V.Ivanchenko - first implementation on base of 
//                         J.P Wellisch class G4PiNuclearCrossSection
// 22.01.07 V.Ivanchenko - add cross section interfaces with Z and A
// 05.03.07 V.Ivanchenko - fix weight for interpolation
// 13.03.07 V.Ivanchenko - cleanup at low energies
// 11.09.09 V.Ivanchenko - fixed bug in interpolation
//

#include "G4UPiNuclearCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PhysicsTable.hh"
#include "G4NistManager.hh"
#include "G4HadronicException.hh"

G4UPiNuclearCrossSection::G4UPiNuclearCrossSection()
 : G4VCrossSectionDataSet("G4UPiNuclearCrossSection")
{
  isInitialized = false;
  piPlusElastic = piPlusInelastic = piMinusElastic = piMinusInelastic = 0;
  piPlus  = G4PionPlus::PionPlus();
  piMinus = G4PionMinus::PionMinus();

  NZ = 16;
  aPower  = 0.75;
  elow    = 20.0*MeV;
  elowest = MeV;
  G4NistManager* nist = G4NistManager::Instance();
  for(G4int i=1; i<93; ++i) {
    APower[i] = G4Pow::GetInstance()->powA(nist->GetAtomicMassAmu(i),aPower);
  }
}

G4UPiNuclearCrossSection::~G4UPiNuclearCrossSection()
{
  piPlusElastic->clearAndDestroy();
  piPlusInelastic->clearAndDestroy();
  piMinusElastic->clearAndDestroy();
  piMinusInelastic->clearAndDestroy();
  delete piPlusElastic;
  delete piPlusInelastic;
  delete piMinusElastic;
  delete piMinusInelastic;
}

G4bool 
G4UPiNuclearCrossSection::IsElementApplicable(const G4DynamicParticle*, 
					      G4int Z, const G4Material*)
{
  return (1 < Z);
}

G4double
G4UPiNuclearCrossSection::GetElasticCrossSection(const G4DynamicParticle* dp,
                                                 G4int Z, G4int A)
{
  G4double cross = 0.0;
  G4PhysicsTable* table = 0;
  const G4ParticleDefinition* part = dp->GetDefinition();
  if(part == piPlus) { table = piPlusElastic; }
  else if(part == piMinus) { table = piMinusElastic; }
  if(table) {
    cross = Interpolate(Z, A, dp->GetKineticEnergy(),table);
  }
  return cross;
}

G4double
G4UPiNuclearCrossSection::GetInelasticCrossSection(const G4DynamicParticle* dp,
                                                   G4int Z, G4int A)
{
  G4double cross = 0.0;
  G4double fact  = 1.0;
  G4double ekin  = dp->GetKineticEnergy();
  G4PhysicsTable* table = 0;
  const G4ParticleDefinition* part = dp->GetDefinition();

  // Coulomb barrier
  if(part == piPlus) {
    if(ekin > elowest) {
      table = piPlusInelastic;
      if(ekin < elow) {
        fact = std::sqrt((ekin-elowest)/(elow-elowest));
        ekin = elow;
      }
    }
  } else if(part == piMinus) {
    table = piMinusInelastic;
    if(ekin < elow) { ekin = elow; }
  }
  if(table) {
    cross = fact*Interpolate(Z, A, ekin, table);
  }
  return cross;
}

G4double G4UPiNuclearCrossSection::Interpolate(
	 G4int Z, G4int A, G4double ekin, G4PhysicsTable* table)
{
  G4double res = 0.0;
  G4int idx;
  G4int iz = Z;
  if(iz > 92) iz = 92;
  for(idx=0; idx<NZ; idx++) {if(theZ[idx] >= iz) break;}
  if(idx >= NZ) idx = NZ - 1;
  G4int iz2 = theZ[idx];
  //  G4cout << "U: iz= " << iz << " iz2= " << iz2 << "  " 
  //  << APower[iz] << "  " << APower[iz2]<<G4endl;
  G4double x2 = (((*table)[idx])->Value(ekin))*APower[iz]/APower[iz2];

  // use only one Z
  if(iz >= theZ[idx] || idx == 0) {
    res = x2;

  // Interpolation between Z
  } else {

    G4int iz1 = theZ[idx-1];
    G4double x1 = (((*table)[idx-1])->Value(ekin))*APower[iz]/APower[iz1];
    G4double w1 = G4double(A) - theA[idx-1];
    G4double w2 = theA[idx] - G4double(A);
    res = (w1*x2 + w2*x1)/(w1 + w2); 
  }
  return res;
}

void G4UPiNuclearCrossSection::AddDataSet(const G4String& p, 
					  const G4double* tot, 
					  const G4double* in, 
					  const G4double* e, 
					  G4int n)
{
  G4LPhysicsFreeVector* pvin = new G4LPhysicsFreeVector(n,e[0]*GeV,e[n-1]*GeV);
  //pvin->SetSpline(true);
  G4LPhysicsFreeVector* pvel = new G4LPhysicsFreeVector(n,e[0]*GeV,e[n-1]*GeV);
  //pvel->SetSpline(true);
  for(G4int i=0; i<n; ++i) { 
    pvin->PutValues(i,e[i]*GeV,in[i]*millibarn); 
    pvel->PutValues(i,e[i]*GeV,std::max(0.0,(tot[i]-in[i])*millibarn)); 
  }
  if(p == "pi+") {
    piPlusInelastic->push_back(pvin);
    piPlusElastic->push_back(pvel);
  } else {
    piMinusInelastic->push_back(pvin);
    piMinusElastic->push_back(pvel);
  }
} 

void G4UPiNuclearCrossSection::DumpPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == piPlus) {
    G4cout << "### G4UPiNuclearCrossSection Elastic data for pi+" << G4endl;
    G4cout << *piPlusElastic << G4endl;
    G4cout << "### G4UPiNuclearCrossSection Inelastic data for pi+" << G4endl;
    G4cout << *piPlusInelastic << G4endl;
  } else if(&p == piMinus) {
    G4cout << "### G4UPiNuclearCrossSection Elastic data for pi-" << G4endl;
    G4cout << *piMinusElastic << G4endl;
    G4cout << "### G4UPiNuclearCrossSection Inelastic data for pi-" << G4endl;
    G4cout << *piMinusInelastic << G4endl;
  }
}

void G4UPiNuclearCrossSection::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(isInitialized) { return; }
  if(&p != piPlus && &p != piMinus) { 
    throw G4HadronicException(__FILE__, __LINE__,"Is applicable only for pions");
    return;
  }
  isInitialized = true;

  const G4int n = 16;
  const G4int iz[n] = {2,4,6,7,8,11,13,20,26,29,42,48,50,74,82,92};
  NZ = n;
  theZ.reserve(n);
  theA.reserve(n);

  G4NistManager* nist = G4NistManager::Instance();
  G4int i;
  for(i=0; i<n; ++i) {
    theZ.push_back(iz[i]);
    theA.push_back(nist->GetAtomicMassAmu(iz[i]));
  }

  piPlusElastic    = new G4PhysicsTable();
  piPlusInelastic  = new G4PhysicsTable();
  piMinusElastic   = new G4PhysicsTable();
  piMinusInelastic = new G4PhysicsTable();

  const G4double e1[38] = {
    0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15, 0.16, 
    0.17, 0.18, 0.19, 0.2,  0.22,0.24, 0.26, 0.28, 0.3,  0.35, 
    0.4,  0.45, 0.5,  0.55, 0.6, 0.7,  0.8,  0.9,  1,    2, 
    3,    5,    10,   20,   50,  100,  500, 1000};
  const G4double e2[39] = {
    0.02, 0.04, 0.06, 0.08, 0.1,  0.11, 0.12, 0.13, 0.14, 0.15, 
    0.16, 0.17, 0.18, 0.2,  0.22, 0.24, 0.26, 0.28, 0.3,  0.35, 
    0.4,  0.45, 0.5,  0.55, 0.575,0.6,  0.7,  0.8,  0.9,  1, 
    2,    3,    5,   10,   20,    50,   100,  500,  1000};
  const G4double e3[31] = {
    0.02, 0.04, 0.06, 0.08, 0.1,  0.12, 0.14, 0.16, 0.18, 0.2, 
    0.22, 0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.6,  0.7,  0.8, 
    0.9,  1,    2,    3,    5,    10,    20,   50,   100, 500, 1000};
  const G4double e4[32] = {
    0.02, 0.04, 0.06, 0.08, 0.1,  0.12, 0.14, 0.16, 0.18,  0.2, 
    0.22, 0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55,  0.6,   0.7,  
    0.8,  0.9,    1,    2,    3,    5,   10,   20,   50,  100,  500,  1000};
  const G4double e5[34] = {
   0.02,  0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.12, 0.14, 
   0.16, 0.18,  0.2,  0.22, 0.25, 0.3,  0.35,  0.4, 0.45,  0.5,  
   0.6,  0.7,   0.8,  0.9,    1,    2,    3,   5,   10,    20, 50, 100, 500, 1000};
  const G4double e6[35] = {
   0.02,  0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.12, 0.14, 
   0.16, 0.18,  0.2,  0.22, 0.25,  0.3, 0.35,  0.4, 0.45,  0.5, 
   0.55,  0.6,  0.7,  0.8,  0.9,    1,    2,    3,    5,   10, 20, 50, 100, 500, 1000};

  const G4double he_t[38] = {
   40,  70, 108,   152, 208,   276, 300, 320,   329, 333, 
  332, 328, 322,   310, 288,   260, 240, 216,   196, 144, 
  125, 112, 108.5, 109, 110.5, 117, 123, 128.5, 135, 110, 
   96,  87,  85,  83.5, 83.5, 83.5, 83.5, 83.5};
  const G4double he_in[38] = {
   18,  38,  62,    98, 136,   176, 190, 200,   209, 212, 
  212, 208, 204,   196, 176,   164, 150, 134,   124, 97.5, 
   90,  85,82.5,  83.5, 86.5,   93, 97.5,100,   102, 83, 
   77,  75,  74,  72.5, 72.5, 72.5, 72.5, 72.5};
  const G4double be_m_t[38] = {
  150, 210, 294,   396, 520,   600, 623, 635,   642, 640, 
  630, 615, 600,   576, 540,   504, 470, 435,   400, 340, 
  294, 258, 236,   230, 233,   244, 257, 270,   276, 250, 
  230, 215, 205,   194, 188,   186, 186, 186};
  const G4double be_m_in[38] = {
   90, 126, 177,   240, 320,   380, 400, 410,   414, 410, 
  400, 387, 371,   360, 333,   312, 285, 260,   237, 216, 
  198, 187, 182,   180, 182,   187, 193, 203,   207, 179, 
  172, 165, 159,   155, 144,   144, 144, 144};
  const G4double be_p_t[38] = {
   96, 150, 222,   320, 430,   514, 545, 565,   574, 574, 
  564, 552, 535,   522, 490,   462, 432, 398,   367, 314, 
   276, 248, 232,  230, 233,   244, 257, 270,   276, 250, 
  230, 215, 205,   194, 188,   186, 186, 186};
  const G4double be_p_in[38] = {
   60,  95, 142,   194, 262,   319, 345, 361,   364, 364, 
  354, 350, 330,   319, 298,   280, 258, 237,   216, 200, 
  189, 183, 182,   180, 182,   187, 193, 203,   207, 179, 
  172, 165, 159,   155, 144,   144, 144, 144};

  const G4double c_m_t[39] = {
  204, 260, 366,   517, 630,   673, 694, 704, 710, 711, 
  706, 694, 676,   648, 616,   584, 548, 518, 489, 426, 
  376, 342, 323,   310, 312,   313, 319, 333, 342, 348, 
  310, 290, 268,   250, 245,   237, 234, 234, 234};
  const G4double c_m_in[39] = {
  128, 160, 224,   315, 388,   416, 430, 438, 444, 445, 
  440, 432, 416,   400, 380,   354, 320, 304, 288, 264, 
  246, 240, 233,   232, 233,   234, 238, 246, 252, 256, 
  220, 210, 198,   187, 183,   176, 174, 174, 174};
  const G4double c_p_t[39] = {
  140, 192, 294,   428, 594,   642, 662, 687, 685, 688, 
  684, 672, 656,   630, 598,   567, 533, 504, 474, 416, 
  369, 336, 319,   310, 312,   313, 319, 333, 342, 348, 
  310, 290, 268,   250, 245,   237, 234, 234, 234};
  const G4double c_p_in[39] = {
   94, 132, 184,   260, 370,   398, 408, 420, 426, 428, 
  424, 416, 400,   386, 366,   340, 308, 294, 280, 257, 
  241, 236, 231,   232, 233,   234, 238, 246, 252, 256, 
  220, 210, 198,   187, 183,   176, 174, 174, 174};
  const G4double n_m_t[39] = {
  246, 308, 424,   590, 729,   776, 800, 821, 822, 817, 
  800, 778, 768,   728, 690,   654, 615, 584, 556, 480, 
  430, 393, 373,   367, 368,   370, 375, 388, 390, 397, 
  364, 337, 310,   291, 275,   268, 268, 268, 268};
  const G4double n_m_in[39] = {
  155, 188, 256,   360, 456,   492, 512, 526, 526, 520, 
  504, 491, 475,   450, 425,   396, 376, 360, 340, 300, 
  282, 270, 265,   265, 266,   268, 273, 280, 288, 288, 
  256, 237, 226,   218, 208,   202, 202, 202, 202};
  const G4double n_p_t[39] = {
  150, 212, 328,   500, 680,   735, 762, 781, 782, 779, 
  770, 748, 740,   706, 672,   633, 600, 569, 541, 467, 
  419, 385, 368,   364, 366,   368, 375, 388, 390, 397, 
  364, 337, 310,   291, 275,   268, 268, 268, 268};
  const G4double n_p_in[39] = {
   90, 140, 208,   300, 426,   467, 490, 504, 504, 500, 
  484, 474, 460,   437, 413,   381, 365, 350, 330, 292, 
  276, 267, 263,   264, 265,   267, 273, 280, 288, 288, 
  256, 237, 226,   218, 208,   202, 202, 202, 202};

  const G4double o_m_t[31] = {
  280, 360, 500,   685, 812,   861, 870, 865, 835, 800, 
  755, 700, 600,   537, 493,   468, 441, 436, 443, 449, 
  460, 463, 432,   385, 350,   325, 312, 307, 303, 303, 303};
  const G4double o_m_in[31] = {
  190, 207, 300,   420, 500,   540, 550, 542, 520, 490, 
  460, 423, 360,   339, 321,   314, 312, 314, 319, 324, 
  328, 330, 300,   275, 250,   240, 229, 225, 222, 222, 222};
  const G4double o_p_t[31] = {
  170, 240, 390,   570, 740,   818, 830, 822, 800, 765, 
  725, 675, 585,   525, 483,   458, 444, 447, 453, 449,
  460, 463, 432,   385, 350,   325, 312, 307, 303, 303, 303};
  const G4double o_p_in[31] = {
  100, 145, 240,   340, 470,   518, 530, 522, 505, 477, 
  448, 412, 350,   330, 316,   310, 308, 311, 317, 324,
  328, 330, 300,   275, 250,   240, 229, 225, 222, 222, 222};
  const G4double na_m_t[31] = {
  450, 545, 705,   910, 1020, 1075, 1087, 1080, 1042, 987, 
  943, 885, 790,   700,  650,  610,  585,  575,  585, 595, 
  600, 610, 556,   524,  494,  458,  445,  429,  427, 427, 427};
  const G4double na_m_in[31] = {
  275, 315, 413,   545,  620,  660,  670,  662,  630, 593, 
  570, 520, 465,   420,  410,  395,  390,  400,  410, 418, 
  420, 422, 372,   348,  330,  320,  310,  294,  292, 292, 292};
  const G4double na_p_t[31] = {
  210, 320, 530,   795,  960, 1035, 1050, 1040, 1007, 957, 
  918, 865, 773,   685,  636,  598,  575,  565,  578, 590, 
  598, 610, 556,   524,  494,  458,  445,  429,  427, 427, 427};
  const G4double na_p_in[31] = {
  115, 210, 340,   495,  585,  630,  645,  637,  605, 572, 
  550, 505, 455,   410,  401,  388,  383,  393,  405, 414, 
  418, 422, 372,   348,  330,  320,  310,  294,  292, 292, 292};
  const G4double al_m_t[31] = {
  532, 637, 832,  1057, 1207, 1230, 1210, 1174, 1133, 1095, 
 1038, 970, 890,   807,  750,  710,  675,  665,  670,  673, 
  678, 682, 618,   574,  546,  520,  507,  495,  488,  488, 488};
  const G4double al_m_in[31] = {
  300, 360, 495,   665,  750,  765,  750,  730,  700,  660, 
  615, 570, 520,   490,  470,  450,  448,  450,  450,  452, 
  456, 460, 408,   392,  376,  356,  347,  338,  332,  332, 332};
  const G4double al_p_t[31] = {
  225, 350, 616,   945, 1122, 1175, 1157, 1128, 1088, 1045, 
  988, 935, 870,   787,  730,  690,  660,  652,  660,  668, 
  678, 682, 618,   574,  546,  520,  507,  495,  488,  488, 488};
  const G4double al_p_in[31] = {
  120, 238, 390,   610,  712,  735,  720,  703,  655,  635, 
  590, 550, 505,   475,  455,  438,  440,  445,  445,  450, 
  456, 460, 408,   392,  376,  356,  347,  338,  332,  332, 332};

  const G4double ca_m_t[31] = {
  800, 980, 1240, 1460, 1570, 1600, 1580, 1535, 1475, 1425, 
 1375,1295, 1200, 1083, 1000,  948,  915,  895,  900,  908, 
  915, 922, 856,   795,  740,  705,  682,  660,  660,  660, 660};
  const G4double ca_m_in[31] = {
  470, 550, 620,   860,  955,  960,  920,  860,  820,  780, 
  740, 665, 637,   615,  600,  590,  580,  580,  600,  608,  
  610, 615, 550,   525,  510,  488,  470,  450,  450,  450,  450};
  const G4double ca_p_t[31] = {
  275, 445, 790,  1195, 1440, 1485, 1475, 1435, 1385, 1335, 
 1295,1245,1160,  1050,  970,  923,  895,  877,  887,  897, 
  904, 913, 855,   795,  740,  705,  682,  660,  660,  660, 660};
  const G4double ca_p_in[31] = {
  160, 315, 500,   745,  870,  905,  900,  860,  810,  770, 
  740, 710, 640,   617,  595,  585,  575,  575,  590,  600, 
  602, 608, 550,   525,  510,  488,  470,  450,  450,  450,  450};
  // ca data may have typo

  const G4double fe_m_t[32] = {
 1175, 1363, 1670, 1950, 2050, 2040, 1975, 1886, 1834, 1773, 
 1720, 1635, 1474, 1380, 1269, 1225, 1182, 1162, 1159, 1162, 
 1178, 1190, 1197, 1102, 1135,  975,  945,  925,  905,  905,  
  905,  905};
  const G4double fe_m_in[32] = {
  625,  725,  910, 1180, 1275, 1250, 1200, 1150, 1100, 1040,  
  995,  925,  825,  810,  780,  760,  745,  740,  740,  740,  
  750,  760,  765,  690,  660,  635,  615,  600,  585,  585,  
  585,  585};
  const G4double fe_p_t[32] =  {
  330,  575, 1010, 1500, 1837, 1875, 1820, 1751, 1691, 1636, 
 1690, 1450, 1396, 1305, 1219, 1190, 1148, 1138, 1134, 1144, 
  1163, 1175, 1183, 1198, 1135, 975,  945,  925,  905,  905,  
  905,  905};
  const G4double fe_p_in[32] = {
  210,  410,  707, 1010, 1125, 1150, 1100, 1070, 1010,  960,  
  920,  776,  780,  760,  750,  740,  720,  725,  725,  730,  
  740,  750,  755,  690,  660,  635,  615,  600,  585,  585,  
  585,  585};
  const G4double cu_m_t[32] = {
 1400, 1600, 1875, 2088, 2200, 2220, 2175, 2125, 2075, 2012, 
 1950, 1855, 1670, 1530, 1430, 1370, 1315, 1315, 1315, 1330, 
 1345, 1360, 1365, 1250, 1185, 1128, 1070, 1035, 1010, 1010, 
 1010, 1010};
  const G4double cu_m_in[32] = {
  725,  840, 1020, 1200, 1295, 1300, 1267, 1240, 1213, 1175, 
 1125, 1042,  950,  900,  860,  840,  830,  832,  835,  840,  
  850,  860,  865,  785,  735,  705,  680,  650,  630,  630,  
  630,  630};
  const G4double cu_p_t[32] = {
  355,  605, 1120, 1630, 1940, 2010, 2010, 1980, 1925, 1895, 
 1830, 1730, 1585, 1490, 1400, 1340, 1290, 1290, 1290, 1310, 
 1330, 1345, 1350, 1240, 1185, 1128, 1070, 1035, 1010, 1010, 
 1010, 1010};
  const G4double cu_p_in[32] = {
  230,  425,  780,  1025, 1155, 1190, 1190, 1180, 1125, 1100, 
 1050, 1000,  900,   870,  835,  815,  810,  812,  815,  825,  
  840,  850,  855,   780,  735,  705,  680,  650,  630,  630,  
  630,  630};

  const G4double mo_m_t[34] = {
  2430, 2610, 2710, 2790, 2880, 2940, 2965, 2970, 2970, 2920, 
  2840, 2720, 2570, 2500, 2365, 2200, 2050, 1926, 1825, 1768, 
  1749, 1750, 1778, 1789, 1808, 1690, 1645, 1530, 1492, 1450, 
  1425, 1425, 1425, 1425};
  const G4double mo_m_in[34] = {
   925, 1125, 1250, 1375, 1500, 1600, 1680, 1750, 1770, 1730, 
  1660, 1580, 1500, 1450, 1330, 1250, 1190, 1140, 1100, 1075, 
  1075, 1070, 1088, 1095, 1110, 1035, 1005,  940,  917,  880, 
   860,  860,  860, 860};
  const G4double mo_p_t[34] = {
   410,  730, 1110, 1530, 1920, 2200, 2385, 2520, 2600, 2630, 
  2575, 2470, 2320, 2285, 2185, 2053, 1945, 1852, 1776, 1719, 
  1710, 1716, 1746, 1759, 1778, 1675, 1645, 1530, 1492, 1450, 
  1425, 1425, 1425, 1425};
  const G4double mo_p_in[34] = {
   270,  540,  825,  975, 1140, 1285, 1400, 1480, 1555, 1580, 
  1525, 1470, 1360, 1340, 1255, 1160, 1120, 1085, 1060, 1045, 
  1045, 1045, 1065, 1075, 1090, 1025, 1005,  940,  917,  880, 
   860,  860,  860, 860};
  const G4double cd_m_t[34] = {
  3060, 3125, 3170, 3220, 3255, 3280, 3290, 3260, 3270, 3200, 
  3120, 3080, 3090, 2920, 2810, 2640, 2362, 2230, 2115, 2050, 
  2020, 2025, 2040, 2070, 2100, 1900, 1795, 1740, 1675, 1645, 
  1625, 1620, 1620, 1620};
  const G4double cd_m_in[34]= {
  1025, 1275, 1440, 1625, 1740, 1800, 1880, 1920, 1980, 1920, 
  1850, 1810, 1720, 1650, 1560, 1450, 1330, 1290, 1245, 1210, 
  1200, 1200, 1205, 1205, 1230, 1130, 1085, 1060, 1000,  985, 
   975,  970,  970, 970};
  const G4double cd_p_t[34] = {
   455,  780, 1170, 1700, 2120, 2400, 2600, 2720, 2820, 2840, 
  2800, 2760, 2720, 2640, 2560, 2450, 2252, 2130, 2035, 1985, 
  1970, 1975, 2005, 2035, 2070, 1880, 1795, 1740, 1675, 1645, 
  1625, 1620, 1620, 1620};
  const G4double cd_p_in[34] = {
   310,  580,  880, 1060, 1270, 1400, 1530, 1610, 1660, 1680, 
  1640, 1600, 1560, 1500, 1430, 1330, 1280, 1230, 1200, 1180, 
  1170, 1175, 1180, 1180, 1210, 1120, 1085, 1060, 1000,  985, 
   975,  970,  970, 970};

  const G4double sn_m_t[35] = {
  3000, 3180, 3250, 3300, 3300, 3410, 3470, 3450, 3410, 3350, 
  3280, 3200, 3120, 3050, 2900, 2630, 2500, 2325, 2190, 2100, 
  2060, 2055, 2055, 2055, 2067, 2085, 2000, 1900, 1835, 1770, 
  1720, 1700, 1695, 1695, 1695};
  const G4double sn_m_in[35] = {
  1050, 1350, 1520, 1650, 1800, 1980, 2070, 2120, 2090, 2050, 
  1980, 1920, 1830, 1770, 1670, 1500, 1435, 1350, 1300, 1230, 
  1220, 1235, 1235, 1235, 1237, 1240, 1160, 1120, 1090, 1065, 
  1040, 1020, 1015, 1015, 1015};
  const G4double sn_p_t[35] =  { 
   465,  800, 1200, 1760, 2170, 2480, 2730, 2885, 2970, 2980, 
  2970, 2890, 2840, 2790, 2620, 2450, 2335, 2205, 2080, 2020, 
  2010, 1990, 1990, 2015, 2030, 2045, 1980, 1890, 1835, 1770, 
  1720, 1700, 1695, 1695, 1695};
  const G4double sn_p_in[35] = { 
   315,  590,  880, 1220, 1460, 1580, 1700, 1770, 1810, 1810, 
  1800, 1730, 1680, 1630, 1530, 1400, 1335, 1270, 1210, 1180, 
  1190, 1190, 1190, 1205, 1210, 1210, 1150, 1115, 1090, 1065, 
  1040, 1020, 1015, 1015, 1015};
  const G4double w_m_t[35] = {
  5200, 5115, 5025, 4975, 4900, 4850, 4780, 4725, 4600, 4490, 
  4355, 4255, 4125, 4040, 3830, 3580, 3330, 3110, 2955, 2860, 
  2852, 2845, 2885, 2900, 2915, 2940, 2800, 2660, 2570, 2490, 
  2460, 2425, 2420, 2420, 2420};
  const G4double w_m_in[35] = {
  1450, 1850, 2100, 2350, 2550, 2700, 2825, 2900, 2850, 2750, 
  2630, 2525, 2400, 2300, 2200, 2070, 1880, 1770, 1715, 1680, 
  1680, 1680, 1685, 1690, 1700, 1720, 1635, 1560, 1530, 1460, 
  1440, 1410, 1410, 1410, 1410};
  const G4double w_p_t[35] = { 
   480,  900, 1500, 2350, 3020, 3420, 3650, 3775, 3875, 3830, 
  3750, 3700, 3630, 3550, 3550, 3290, 3070, 2890, 2840, 2730, 
  2725, 2720, 2770, 2805, 2828, 2865, 2770, 2640, 2570, 2490,
  2460, 2425, 2420, 2420, 2420};
  const G4double w_p_in[35] = { 
   325,  680,  990, 1500, 1850, 2150, 2250, 2300, 2350, 2330, 
  2280, 2230, 2200, 2120, 2130, 1900, 1780, 1670, 1635, 1600, 
  1602, 1605, 1610, 1615, 1630, 1660, 1620, 1550, 1530, 1460,
  1440, 1410, 1410, 1410, 1410};

  const G4double pb_m_t[35] = {
  5890, 5700, 5610, 5580, 5550, 5480, 5400, 5300, 5100, 4930, 
  4750, 4600, 4400, 4280, 4170, 3915, 3650, 3470, 3260, 3150, 
  3120, 3070, 3085, 3100, 3120, 3160, 3070, 2930, 2820, 2750, 
  2710, 2655, 2640, 2640, 2640};
  const G4double pb_m_in[35] = {
  1575, 2025, 2300, 2575, 2850, 3000, 3115, 3180, 3080, 2940, 
  2800, 2670, 2550, 2450, 2370, 2220, 2110, 2000, 1920, 1880, 
  1850, 1800, 1805, 1810, 1820, 1840, 1800, 1720, 1640, 1620, 
  1570, 1530, 1530, 1530, 1530};
  const G4double pb_p_t[35] = { 
   515,  940, 1500, 2400, 3270, 3750, 4050, 4140, 4260, 4200, 
  4080, 3990, 3990, 3810, 3730, 3520, 3370, 3186, 3110, 3010, 
  2990, 2985, 3005, 3020, 3040, 3080, 3020, 2905, 2790, 2750,
  2710, 2655, 2640, 2640, 2640};
  const G4double pb_p_in[35] = { 
   348,  707, 1040, 1650, 2100, 2400, 2580, 2640, 2650, 2520, 
  2410, 2300, 2250, 2190, 2130, 2000, 1930, 1870, 1830, 1790, 
  1770, 1765, 1775, 1780, 1790, 1800, 1775, 1710, 1620, 1620,
  1570, 1530, 1530, 1530, 1530};
  const G4double u_m_t[35] = {
  7080, 6830, 6650, 6530, 6400, 6280, 6100, 5840, 5660, 5520, 
  5330, 5160, 4990, 4810, 4630, 4323, 4130, 3870, 3700, 3550, 
  3490, 3465, 3467, 3475, 3495, 3515, 3440, 3360, 3150, 3040, 
  2985, 2955, 2940, 2940, 2940};
  const G4double u_m_in[35] = {
  1740, 2220, 2500, 2820, 3080, 3300, 3420, 3500, 3420, 3330, 
  3200, 3060, 2940, 2850, 2710, 2470, 2380, 2250, 2160, 2080, 
  2040, 2045, 2047, 2050, 2055, 2060, 2010, 1980, 1830, 1780, 
  1735, 1710, 1700, 1700, 1700};
  const G4double u_p_t[35] = { 
   485,  960, 1580, 2700, 3550, 4050, 4320, 4420, 4620, 4660, 
  4580, 4470, 4350, 4295, 4187, 3938, 3755, 3573, 3450, 3342, 
  3310, 3295, 3310, 3330, 3375, 3405, 3350, 3338, 3135, 3040,
  2985, 2955, 2940, 2940, 2940};
  const G4double u_p_in[35] = { 
   334,  720, 1020, 1560, 2100, 2300, 2550, 2700, 2880, 2880, 
  2760, 2660, 2550, 2510, 2430, 2270, 2130, 2060, 2000, 1970, 
  1950, 1950, 1960, 1960, 1970, 1980, 1950, 1978, 1830, 1780,
  1735, 1710, 1700, 1700, 1700};

  AddDataSet("pi-",he_t,   he_in,  e1, 38);
  AddDataSet("pi+",he_t,   he_in,  e1, 38);
  AddDataSet("pi-",be_m_t, be_m_in, e1, 38);
  AddDataSet("pi+",be_p_t, be_p_in, e1, 38);
  AddDataSet("pi-",c_m_t,  c_m_in,  e2, 39);
  AddDataSet("pi+",c_p_t,  c_p_in,  e2, 39);
  AddDataSet("pi-",n_m_t,  n_m_in,  e2, 39);
  AddDataSet("pi+",n_p_t,  n_p_in,  e2, 39);
  AddDataSet("pi-",o_m_t,  o_m_in,  e3, 31);
  AddDataSet("pi+",o_p_t,  o_p_in,  e3, 31);
  AddDataSet("pi-",na_m_t, na_m_in, e3, 31);
  AddDataSet("pi+",na_p_t, na_p_in, e3, 31);
  AddDataSet("pi-",al_m_t, al_m_in, e3, 31);
  AddDataSet("pi+",al_p_t, al_p_in, e3, 31);
  AddDataSet("pi-",ca_m_t, ca_m_in, e3, 31);
  AddDataSet("pi+",ca_p_t, ca_p_in, e3, 31);
  AddDataSet("pi-",fe_m_t, fe_m_in, e4, 32);
  AddDataSet("pi+",fe_p_t, fe_p_in, e4, 32);
  AddDataSet("pi-",cu_m_t, cu_m_in, e4, 32);
  AddDataSet("pi+",cu_p_t, cu_p_in, e4, 32);
  AddDataSet("pi-",mo_m_t, mo_m_in, e5, 34);
  AddDataSet("pi+",mo_p_t, mo_p_in, e5, 34);
  AddDataSet("pi-",cd_m_t, cd_m_in, e5, 34);
  AddDataSet("pi+",cd_p_t, cd_p_in, e5, 34);
  AddDataSet("pi-",sn_m_t, sn_m_in, e6, 35);
  AddDataSet("pi+",sn_p_t, sn_p_in, e6, 35);
  AddDataSet("pi-",w_m_t,  w_m_in,  e6, 35);
  AddDataSet("pi+",w_p_t,  w_p_in,  e6, 35);
  AddDataSet("pi-",pb_m_t, pb_m_in, e6, 35);
  AddDataSet("pi+",pb_p_t, pb_p_in, e6, 35);
  AddDataSet("pi-",u_m_t,  u_m_in,  e6, 35);
  AddDataSet("pi+",u_p_t,  u_p_in,  e6, 35);
}

void G4UPiNuclearCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4UPiNuclearCrossSection calculates the total, elastic and\n"
          << "inelastic cross sections for pion scattering from nuclei\n"
          << "heavier than hydrogen.  It is based on the Barashenkov\n"
          << "parameterization and is valid for all incident energies.\n";
}

