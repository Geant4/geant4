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
//
// $Id: Tst50ElectronEEDLrange.hh,v 1.2 2003-05-17 18:11:52 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG         Designed for modular Physics List with
// CSDA and StoppingPower test conditions
//
// -------------------------------------------------------------------


#ifndef TST50ELECTRONEEDLrange_HH
#define TST50ELECTRONEEDLrange_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ElectronEEDLrange : public G4VPhysicsConstructor {

public: 

  Tst50ElectronEEDLrange(const G4String& name = "electron-eedl-range");
  
  virtual ~Tst50ElectronEEDLrange();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif








