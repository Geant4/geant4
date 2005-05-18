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
// $Id: HadrontherapyIonLowE.hh,v 1.1 2005-05-18 08:00:00 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it

#ifndef HADRONTHERAPYIONLOWE_HH
#define HADRONTHERAPYIONLOWE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class HadrontherapyIonLowE: public G4VPhysicsConstructor {

public: 

  HadrontherapyIonLowE(const G4String& name = "ion-LowE");
  
  virtual ~HadrontherapyIonLowE();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  virtual void ConstructProcess();
};
#endif

