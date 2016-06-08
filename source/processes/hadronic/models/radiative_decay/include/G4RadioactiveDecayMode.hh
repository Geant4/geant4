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
#ifndef G4RadioactiveDecayMode_h
#define G4RadioactiveDecayMode_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecayMode.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// Identifies a type G4RadioactiveDecayMode to assign a specific decay mode
// description to decay channels.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
////////////////////////////////////////////////////////////////////////////////
//
enum G4RadioactiveDecayMode
  {ERROR=-1, IT=0, BetaMinus=1, BetaPlus=2, KshellEC=3, LshellEC=4, MshellEC=5, Alpha=6};

G4std::istream &operator >> (G4std::istream &s, G4RadioactiveDecayMode &q);

////////////////////////////////////////////////////////////////////////////////
#endif







