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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4GlaubAADataSet.cc
/// \brief Implementation of the G4GlaubAADataSet class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4GlaubAADataSet.cc
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET


#include "G4GlaubAADataSet.hh"

#include "G4DPMJET2_5Interface.hh"

#include <iomanip>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//
// G4GlaubAADataSet
//
// Constructor simply resets all variables to zero.
//
G4GlaubAADataSet::G4GlaubAADataSet() : G4VGlauberDataSet()
{
  ZP    = -1;
  ZT    = -1;

  DefineAZStabilityLine();
}
///////////////////////////////////////////////////////////////////////////////
//
// ~G4GlaubAADataSet
//
// If you thought the contructor was boring, the destructor is even worse!.
// It doesn't do anything.
//
G4GlaubAADataSet::~G4GlaubAADataSet()
{}
////////////////////////////////////////////////////////////////////////////////
//
G4bool G4GlaubAADataSet::CreateGlauberData (const G4int , const G4int )
{
//
//
// Dummy member function.
//
  return false;
}
void G4GlaubAADataSet::DefineAZStabilityLine ()
{
        stabZ[0]  = 0;
        stabZ[1]  = 1;
        stabZ[2]  = 1;
        stabZ[3]  = 2;
        stabZ[4]  = 2;
        stabZ[5]  = 3;
        stabZ[6]  = 3;
        stabZ[7]  = 3;
        stabZ[8]  = 3;
        stabZ[9]  = 4;
        stabZ[10] = 5;
        stabZ[11] = 5;
        stabZ[12] = 6;
        stabZ[13] = 6;
        stabZ[14] = 7;
        stabZ[15] = 7;
        stabZ[16] = 8;
        stabZ[17] = 8;
        stabZ[18] = 8;
        stabZ[19] = 9;
        stabZ[20] = 10;
        stabZ[21] = 10;
        stabZ[22] = 10;
        stabZ[23] = 11;
        stabZ[24] = 12;
        stabZ[25] = 12;
        stabZ[26] = 12;
        stabZ[27] = 13;
        stabZ[28] = 13;
        stabZ[29] = 14;
        stabZ[30] = 14;
        stabZ[31] = 15;
        stabZ[32] = 16;
        stabZ[33] = 16;
        stabZ[34] = 16;
        stabZ[35] = 17;
        stabZ[36] = 17;
        stabZ[37] = 17;
        stabZ[38] = 18;
        stabZ[39] = 19;
        stabZ[40] = 19;
        stabZ[41] = 19;
        stabZ[42] = 20;
        stabZ[43] = 20;
        stabZ[44] = 20;
        stabZ[45] = 21;
        stabZ[46] = 21;
        stabZ[47] = 22;
        stabZ[48] = 21;
        stabZ[49] = 22;
        stabZ[50] = 23;
        stabZ[51] = 23;
        stabZ[52] = 24;
        stabZ[53] = 24;
        stabZ[54] = 25;
        stabZ[55] = 25;
        stabZ[56] = 26;
        stabZ[57] = 26;
        stabZ[58] = 27;
        stabZ[59] = 27;
        stabZ[60] = 28;
        stabZ[61] = 28;
        stabZ[62] = 28;
        stabZ[63] = 29;
        stabZ[64] = 29;
        stabZ[65] = 29;
        stabZ[66] = 30;
        stabZ[67] = 30;
        stabZ[68] = 30;
        stabZ[69] = 31;
        stabZ[70] = 31;
        stabZ[71] = 31;
        stabZ[72] = 32;
        stabZ[73] = 32;
        stabZ[74] = 32;
        stabZ[75] = 33;
        stabZ[76] = 33;
        stabZ[77] = 34;
        stabZ[78] = 35;
        stabZ[79] = 35;
        stabZ[80] = 35;
        stabZ[81] = 35;
        stabZ[82] = 35;
        stabZ[83] = 36;
        stabZ[84] = 37;
        stabZ[85] = 37;
        stabZ[86] = 37;
        stabZ[87] = 37;
        stabZ[88] = 38;
        stabZ[89] = 39;
        stabZ[90] = 40;
        stabZ[91] = 40;
        stabZ[92] = 41;
        stabZ[93] = 41;
        stabZ[94] = 41;
        stabZ[95] = 42;
        stabZ[96] = 42;
        stabZ[97] = 42;
        stabZ[98] = 43;
        stabZ[99] = 44;
        stabZ[100]= 43;
        stabZ[101]= 44;
        stabZ[102]= 45;
        stabZ[103]= 45;
        stabZ[104]= 45;
        stabZ[105]= 46;
        stabZ[106]= 47;
        stabZ[107]= 47;
        stabZ[108]= 47;
        stabZ[109]= 47;
        stabZ[110]= 47;
        stabZ[111]= 48;
        stabZ[112]= 49;
        stabZ[113]= 49;
        stabZ[114]= 49;
        stabZ[115]= 49;
        stabZ[116]= 49;
        stabZ[117]= 50;
        stabZ[118]= 50;
        stabZ[119]= 50;
        stabZ[120]= 51;
        stabZ[121]= 51;
        stabZ[122]= 51;
        stabZ[123]= 51;
        stabZ[124]= 52;
        stabZ[125]= 52;
        stabZ[126]= 53;
        stabZ[127]= 53;
        stabZ[128]= 53;
        stabZ[129]= 54;
        stabZ[130]= 54;
        stabZ[131]= 54;
        stabZ[132]= 55;
        stabZ[133]= 55;
        stabZ[134]= 55;
        stabZ[135]= 56;
        stabZ[136]= 55;
        stabZ[137]= 56;
        stabZ[138]= 57;
        stabZ[139]= 57;
        stabZ[140]= 58;
        stabZ[141]= 59;
        stabZ[142]= 59;
        stabZ[143]= 60;
        stabZ[144]= 61;
        stabZ[145]= 60;
        stabZ[146]= 60;
        stabZ[147]= 62;
        stabZ[148]= 61;
        stabZ[149]= 62;
        stabZ[150]= 61;
        stabZ[151]= 63;
        stabZ[152]= 63;
        stabZ[153]= 63;
        stabZ[154]= 63;
        stabZ[155]= 64;
        stabZ[156]= 64;
        stabZ[157]= 64;
        stabZ[158]= 65;
        stabZ[159]= 65;
        stabZ[160]= 65;
        stabZ[161]= 66;
        stabZ[162]= 67;
        stabZ[163]= 66;
        stabZ[164]= 67;
        stabZ[165]= 67;
        stabZ[166]= 68;
        stabZ[167]= 68;
        stabZ[168]= 69;
        stabZ[169]= 69;
        stabZ[170]= 69;
        stabZ[171]= 70;
        stabZ[172]= 70;
        stabZ[173]= 70;
        stabZ[174]= 70;
        stabZ[175]= 71;
        stabZ[176]= 71;
        stabZ[177]= 72;
        stabZ[178]= 72;
        stabZ[179]= 72;
        stabZ[180]= 73;
        stabZ[181]= 73;
        stabZ[182]= 74;
        stabZ[183]= 74;
        stabZ[184]= 75;
        stabZ[185]= 75;
        stabZ[186]= 75;
        stabZ[187]= 76;
        stabZ[188]= 76;
        stabZ[189]= 76;
        stabZ[190]= 77;
        stabZ[191]= 77;
        stabZ[192]= 77;
        stabZ[193]= 77;
        stabZ[194]= 78;
        stabZ[195]= 78;
        stabZ[196]= 79;
        stabZ[197]= 79;
        stabZ[198]= 79;
        stabZ[199]= 80;
        stabZ[200]= 80;
        stabZ[201]= 80;
        stabZ[202]= 80;
        stabZ[203]= 81;
        stabZ[204]= 81;
        stabZ[205]= 81;
        stabZ[206]= 82;
        stabZ[207]= 82;
        stabZ[208]= 82;
        stabZ[209]= 83;
        stabZ[210]= 83;
}

#endif
