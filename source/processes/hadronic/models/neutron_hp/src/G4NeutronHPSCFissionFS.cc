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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPSCFissionFS.hh"

  void G4NeutronHPSCFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String aString = "/SC/";
    G4NeutronHPFissionBaseFS::Init(A, Z, dirName, aString);
  }
  
  G4DynamicParticleVector * G4NeutronHPSCFissionFS::ApplyYourself(G4int NNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4NeutronHPSCFissionFS::ApplyYourself +"<<G4endl;
    aResult = G4NeutronHPFissionBaseFS::ApplyYourself(NNeutrons);    
    return aResult;
  }
