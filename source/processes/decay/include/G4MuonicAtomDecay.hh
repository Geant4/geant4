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
// $Id: G4MuonicAtomDecay.hh 78547 2014-01-07 09:40:38Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, 30 August 2016 K.Lynch


#ifndef G4MuonicAtomDecay_h
#define G4MuonicAtomDecay_h 1

#include "G4Decay.hh"

class G4MuonicAtomDecay : public G4Decay
{
 // Class Description
  //  This class is a decay process

  public:
    //  Constructors 
    G4MuonicAtomDecay(const G4String& processName ="MuonicAtomDecay");

    virtual G4bool IsApplicable(const G4ParticleDefinition&) override;
    // returns "true" if the decay process can be applied to
    // the particle type. 
};
#endif










