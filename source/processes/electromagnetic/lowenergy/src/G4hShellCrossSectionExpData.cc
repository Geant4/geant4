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
//
// Author: Simona Saliceti (simona.saliceti@ge.infn.it)
//
// History:
// -----------
// 22 Apr 2004  First committed to cvs
//
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionExpData.cc,v 1.4 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4hShellCrossSectionExpData.hh"
#include <fstream>

G4hShellCrossSectionExpData::G4hShellCrossSectionExpData()
{


  FillVectorValues();

  FillParameterMap(); 

}

G4hShellCrossSectionExpData::~G4hShellCrossSectionExpData()
{ 
  std::map< G4int,std::vector<G4double>*,std::less<G4int> >::iterator pos;
  for (pos = parameterMap.begin(); pos != parameterMap.end(); ++pos)
    {
      std::vector<G4double>* dataSet = (*pos).second;
      delete dataSet;
    }
}

inline void G4hShellCrossSectionExpData::InitializeVector(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3) const
{
  vect.push_back(value1);
  vect.push_back(value2);
  vect.push_back(value3);
}

void G4hShellCrossSectionExpData::FillVectorValues()
{
  InitializeVector(parameter6C, 25965819.,    0.1974268719,  -0.9746342);
  InitializeVector(parameter7N, 14254204.,    0.09584816998, -0.8908825);
  InitializeVector(parameter8O, 14257000.,    0.9582500E-01, -0.8909900);
  InitializeVector(parameter9F,  5628098.,    0.2221402E-01, -0.7857742);
  InitializeVector(parameter10Ne, 3866500.,   0.9298651E-02, -0.7569925);
  InitializeVector(parameter11Na, 2816614.,   0.3161790E-02, -0.7429561);
  InitializeVector(parameter12Mg, 2113709.,   0.9267687E-03, -0.7342144);
  InitializeVector(parameter13Al, 1647127.,   0.2344573E-03, -0.7305032);
  InitializeVector(parameter14Si, 1267490.,   0.5787061E-04, -0.7184879);
  InitializeVector(parameter15P,  1035792.,   0.1101495E-04, -0.7199858);
  InitializeVector(parameter16S,   825871.6,  0.2326387E-05, -0.7095044);
  InitializeVector(parameter17Cl,      13.41653,  -14.66877, -0.7017217);
  InitializeVector(parameter18Ar,  555038.1,  0.6543986E-07, -0.6968985);
  InitializeVector(parameter19K,   469063.8,  0.8876799E-08, -0.6940991);
  InitializeVector(parameter20Ca,  405918.4,  0.9881247E-09, -0.6953619);
  InitializeVector(parameter21Sc,  351611.0,  0.1054022E-09, -0.6950487);
  InitializeVector(parameter22Ti,  313223.2,  0.8541349E-11, -0.6984913);
  InitializeVector(parameter23V,   257752.3,  0.1092866E-11, -0.6841001);
  InitializeVector(parameter24Cr,  230064.8,  0.7874358E-13, -0.6862135);
  InitializeVector(parameter25Mn,  151740.8,  0.7988999E-13, -0.6279536);
  InitializeVector(parameter26Fe,      12.21083,  -35.54018, -0.7140945);
  InitializeVector(parameter27Co,      12.14600,  -38.73178, -0.7229655);
  InitializeVector(parameter28Ni,      12.10570,  -42.25900, -0.7354265);
  InitializeVector(parameter29Cu,      11.87911,  -44.37304, -0.7089269);
  InitializeVector(parameter30Zn,      11.86411,  -48.35208, -0.7252861);
  InitializeVector(parameter31Ga,      11.66042,  -50.73798, -0.7029260);
  InitializeVector(parameter32Ge,      11.66824,  -55.26208, -0.7218891);
  InitializeVector(parameter33As,      11.47804,  -57.94218, -0.7016547);
  InitializeVector(parameter34Se,      11.30017,  -60.70304, -0.6836499);
  InitializeVector(parameter35Br,      11.32571,  -65.77174, -0.7042627);
  InitializeVector(parameter36Kr,      10.44025,  -61.74093, -0.5383769);
  InitializeVector(parameter37Rb,      10.27499,  -64.44025, -0.5236434);
  InitializeVector(parameter38Sr,      10.31786,  -69.83844, -0.5475903);
  InitializeVector(parameter39Y,       10.14157,  -72.55550, -0.5302060);
  InitializeVector(parameter40Zr,      10.03310,  -76.52999, -0.5243382);
  InitializeVector(parameter41Nb,      10.03804,  -81.38335, -0.5404782);
  InitializeVector(parameter42Mo,       9.872166, -84.24110, -0.5245091);
  InitializeVector(parameter43Tc,       9.941685, -90.72401, -0.5503877);
  InitializeVector(parameter44Ru,       9.771618, -93.50847, -0.5334178);
  InitializeVector(parameter45Rh,       9.819192, -99.99145, -0.5543751);
  InitializeVector(parameter46Pd,       9.680034, -103.5450, -0.5387755);
  InitializeVector(parameter47Ag,       9.703977, -109.8104, -0.5550457);
  InitializeVector(parameter48Cd,       9.548004, -113.0565, -0.5399088);
  InitializeVector(parameter49In,       9.378392, -115.8466, -0.5226109);
  InitializeVector(parameter50Sn,       9.406608, -122.6137, -0.5386317);
  InitializeVector(parameter51Sb,       9.254029, -125.8533, -0.5238908);
  InitializeVector(parameter52Te,       9.311451, -133.6049, -0.5440237);
  InitializeVector(parameter53I,        9.158569, -136.7913, -0.5291739);
  InitializeVector(parameter54Xe,       9.011420, -140.1908, -0.5150210);
  InitializeVector(parameter55Cs,       9.112158, -149.6160, -0.5411919);
  InitializeVector(parameter56Ba,       7.696734, -124.4350, -0.3181673);
  InitializeVector(parameter57La,       7.728227, -130.9067, -0.3348193);
  InitializeVector(parameter58Ce,       7.595518, -134.1356, -0.3233585);
  InitializeVector(parameter59Pr,       7.451400, -136.9318, -0.3102762);
  InitializeVector(parameter60Nd,       7.493138, -144.0968, -0.3273972);
  InitializeVector(parameter61Pm,       7.358854, -147.1354, -0.3155824);
  InitializeVector(parameter62Sm,       7.231769, -150.3722, -0.3048225);
  InitializeVector(parameter63Eu,       7.060385, -152.3517, -0.2873322);
  InitializeVector(parameter64Gd,       7.179881, -161.9544, -0.3162093);
  InitializeVector(parameter65Tb,       7.048202, -164.9670, -0.3046019);
  InitializeVector(parameter66Dy,       1.962871, -4.620275, 0.5748519E-02);
  InitializeVector(parameter67Ho,       1.992797, -4.863559, 0.5772170E-02);
  InitializeVector(parameter68Er,       2.013621, -5.068997, 0.5749294E-02);
  InitializeVector(parameter69Tm,       2.026489, -5.243042, 0.5671808E-02);
  InitializeVector(parameter70Yb,       1.946252, -5.047087, 0.5011610E-02);
  InitializeVector(parameter71Lu,       1.968246, -5.253083, 0.5007004E-02);
  InitializeVector(parameter72Hf,       1.980888, -5.421401, 0.4949681E-02);
  InitializeVector(parameter73Ta,       1.925043, -5.310529, 0.4504530E-02);
  InitializeVector(parameter74W,        1.939689, -5.484359, 0.4464923E-02);
  InitializeVector(parameter75Re,       1.951746, -5.645718, 0.4415765E-02);
  InitializeVector(parameter76Os,       2.249802, -6.974045, 0.5965646E-02);
  InitializeVector(parameter77Ir,       2.257475, -7.121614, 0.5873831E-02);
  InitializeVector(parameter78Pt,       2.262231, -7.255240, 0.5769194E-02);
  InitializeVector(parameter79Au,       2.188981, -7.080360, 0.5161316E-02);
  InitializeVector(parameter80Hg,       2.199681, -7.233214, 0.5111361E-02);
  InitializeVector(parameter81Tl,       2.211119, -7.389515, 0.5064001E-02);
  InitializeVector(parameter82Pb,       2.204573, -7.468263, 0.4934075E-02);
  InitializeVector(parameter83Bi,       2.148083, -7.348053, 0.4494417E-02);
  InitializeVector(parameter84Po,       2.154512, -7.479838, 0.4433237E-02);
  InitializeVector(parameter85At,       2.164621, -7.626572, 0.4388465E-02);
  InitializeVector(parameter86Rn,       2.288465, -8.212229, 0.5807703E-02);
  InitializeVector(parameter87Fr,       2.294291, -8.336506, 0.5753172E-02);
  InitializeVector(parameter88Ra,       2.318430, -8.525764, 0.5849121E-02);
  InitializeVector(parameter89Ac,       2.306637, -8.578354, 0.5672658E-02);
  InitializeVector(parameter90Th,       2.254484, -8.479531, 0.5163647E-02);
  InitializeVector(parameter91Pa,       2.246310, -8.540535, 0.5044748E-02);
  InitializeVector(parameter92U,        2.276106, -8.748476, 0.5186345E-02);
}

void G4hShellCrossSectionExpData::FillParameterMap()
{
  parameterMap [6] =  &parameter6C ;
  parameterMap [7] =  &parameter7N ;
  parameterMap [8] =  &parameter8O ;
  parameterMap [9] =  &parameter9F ;
  parameterMap [10] = &parameter10Ne;
  parameterMap [11] = &parameter11Na;
  parameterMap [12] = &parameter12Mg;
  parameterMap [13] = &parameter13Al;
  parameterMap [14] = &parameter14Si;
  parameterMap [15] = &parameter15P ;
  parameterMap [16] = &parameter16S ;
  parameterMap [17] = &parameter17Cl;
  parameterMap [18] = &parameter18Ar;
  parameterMap [19] = &parameter19K ;
  parameterMap [20] = &parameter20Ca;
  parameterMap [21] = &parameter21Sc;
  parameterMap [22] = &parameter22Ti;
  parameterMap [23] = &parameter23V ;
  parameterMap [24] = &parameter24Cr;
  parameterMap [25] = &parameter25Mn;
  parameterMap [26] = &parameter26Fe;
  parameterMap [27] = &parameter27Co;
  parameterMap [28] = &parameter28Ni;
  parameterMap [29] = &parameter29Cu;
  parameterMap [30] = &parameter30Zn;
  parameterMap [31] = &parameter31Ga;
  parameterMap [32] = &parameter32Ge;
  parameterMap [33] = &parameter33As;
  parameterMap [34] = &parameter34Se;
  parameterMap [35] = &parameter35Br;
  parameterMap [36] = &parameter36Kr;
  parameterMap [37] = &parameter37Rb;
  parameterMap [38] = &parameter38Sr;
  parameterMap [39] = &parameter39Y ;
  parameterMap [40] = &parameter40Zr;
  parameterMap [41] = &parameter41Nb;
  parameterMap [42] = &parameter42Mo;
  parameterMap [43] = &parameter43Tc;
  parameterMap [44] = &parameter44Ru;
  parameterMap [45] = &parameter45Rh;
  parameterMap [46] = &parameter46Pd;
  parameterMap [47] = &parameter47Ag;
  parameterMap [48] = &parameter48Cd;
  parameterMap [49] = &parameter49In;
  parameterMap [50] = &parameter50Sn;
  parameterMap [51] = &parameter51Sb;
  parameterMap [52] = &parameter52Te;
  parameterMap [53] = &parameter53I ;
  parameterMap [54] = &parameter54Xe;
  parameterMap [55] = &parameter55Cs; 
  parameterMap [56] = &parameter56Ba;
  parameterMap [57] = &parameter57La;
  parameterMap [58] = &parameter58Ce;
  parameterMap [59] = &parameter59Pr;
  parameterMap [60] = &parameter60Nd;
  parameterMap [61] = &parameter61Pm;
  parameterMap [62] = &parameter62Sm;
  parameterMap [63] = &parameter63Eu;
  parameterMap [64] = &parameter64Gd;
  parameterMap [65] = &parameter65Tb;
  parameterMap [66] = &parameter66Dy;
  parameterMap [67] = &parameter67Ho;
  parameterMap [68] = &parameter68Er;
  parameterMap [69] = &parameter69Tm;
  parameterMap [70] = &parameter70Yb;
  parameterMap [71] = &parameter71Lu;
  parameterMap [72] = &parameter72Hf;
  parameterMap [73] = &parameter73Ta;
  parameterMap [74] = &parameter74W ;
  parameterMap [75] = &parameter75Re;
  parameterMap [76] = &parameter76Os;
  parameterMap [77] = &parameter77Ir;
  parameterMap [78] = &parameter78Pt;
  parameterMap [79] = &parameter79Au;
  parameterMap [80] = &parameter80Hg;
  parameterMap [81] = &parameter81Tl;
  parameterMap [82] = &parameter82Pb;
  parameterMap [83] = &parameter83Bi;
  parameterMap [84] = &parameter84Po;
  parameterMap [85] = &parameter85At;
  parameterMap [86] = &parameter86Rn;
  parameterMap [87] = &parameter87Fr;
  parameterMap [88] = &parameter88Ra;
  parameterMap [89] = &parameter89Ac;
  parameterMap [90] = &parameter90Th;
  parameterMap [91] = &parameter91Pa;
  parameterMap [92] = &parameter92U ;
}

std::vector<G4double>* G4hShellCrossSectionExpData::GetParam(G4int Z)
{   
  return parameterMap[Z];
}
