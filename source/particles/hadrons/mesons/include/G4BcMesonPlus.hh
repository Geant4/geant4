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
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
//      Created,             Hisaya Kurashige, 11 Aug. 2011
// **********************************************************************
// ----------------------------------------------------------------

#ifndef G4BcMesonPlus_h
#define G4BcMesonPlus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                        BcMesonPlus                             ###
// ######################################################################

class G4BcMesonPlus : public G4ParticleDefinition
{
 private:
   static G4BcMesonPlus* theInstance;
   G4BcMesonPlus(){}
   ~G4BcMesonPlus(){}

 public:
   static G4BcMesonPlus* Definition();
   static G4BcMesonPlus* BcMesonPlusDefinition();
   static G4BcMesonPlus* BcMesonPlus();
};

#endif
