//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4HadParticleCodes_hh
#define G4HadParticleCodes_hh

enum
{
  NeutronPC = 2112,
  ProtonPC = 2212,
  
  DeltamPC = 1114,
  Delta0PC = 2114,
  DeltapPC = 2214,
  DeltappPC = 2224,
  
  Dm_1600PC = 31114,
  D0_1600PC = 32114,
  Dp_1600PC = 32214,
  Dpp_1600PC = 32224,
  
  Dm_1620PC = 1112,
  D0_1620PC = 1212,
  Dp_1620PC = 2122,
  Dpp_1620PC = 2222,
  
  Dm_1700PC = 11114,
  D0_1700PC = 12114,
  Dp_1700PC = 12214,
  Dpp_1700PC = 12224,
  
  Dm_1900PC = 11112,
  D0_1900PC = 11212,
  Dp_1900PC = 12122,
  Dpp_1900PC = 12222,

  Dm_1905PC = 1116,
  D0_1905PC = 1216,
  Dp_1905PC = 2126,
  Dpp_1905PC = 2226,

  Dm_1910PC = 21112,
  D0_1910PC = 21212,
  Dp_1910PC = 22122,
  Dpp_1910PC = 22222,

  Dm_1920PC = 21114,
  D0_1920PC = 22114,
  Dp_1920PC = 22214,
  Dpp_1920PC = 22224,

  Dm_1930PC = 11116,
  D0_1930PC = 11216,
  Dp_1930PC = 12126,
  Dpp_1930PC = 12226,

  Dm_1950PC = 1118,
  D0_1950PC = 2118,
  Dp_1950PC = 2218,
  Dpp_1950PC = 2228,

  N1400pPC = 12212,
  N1400nPC = 12112,
  
  N1520pPC = 2124,
  N1520nPC = 1214,

  N1535pPC = 22212,
  N1535nPC = 22112,

  N1650pPC = 32212,
  N1650nPC = 32112,

  N1675pPC = 2216,
  N1675nPC = 2116,

  N1680pPC = 12216,
  N1680nPC = 12116,

  N1700pPC = 22124,
  N1700nPC = 21214,

  N1710pPC = 42212,
  N1710nPC = 42112,

  N1720pPC = 32124,
  N1720nPC = 31214,

  N1900pPC = 42124,
  N1900nPC = 41214,

  N1990pPC = 12218,
  N1990nPC = 12118,

  N2090pPC = 52214,
  N2090nPC = 52114,

  N2190pPC = 2128,
  N2190nPC = 1218,

  N2220pPC = 100002210,
  N2220nPC = 100002110,

  N2250pPC = 100012210,
  N2250nPC = 100012110

};

struct D1232
{
  enum
  {
    Dm=DeltamPC,
    D0=Delta0PC,
    Dp=DeltapPC,
    Dpp=DeltappPC
  };
};

#endif
