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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG          Designed for modular Physics List with
// Ziegler model for the Stopping Power and CSDArange and Stopping Power 
//conditions set
//
// -------------------------------------------------------------------
#ifndef TST50PROTONEEDLziegler_HH
#define TST50PROTONEEDLziegler_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ProtonEEDLziegler : public G4VPhysicsConstructor {

public: 

  Tst50ProtonEEDLziegler(const G4String& name = "proton-eedl-ziegler");
  
  virtual ~Tst50ProtonEEDLziegler();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

