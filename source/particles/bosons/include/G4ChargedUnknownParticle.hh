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
// --------------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4ChargedUnknownParticle.hh
//
//      Author:        A.Ribon
// 
//      Creation date: August 2024
//
//      Description:   This class is similar to G4UnknownParticle,
//                     i.e. representing particles with valid PDG code but
//                     unknown to Geant4 - e.g. produced by MC event
//                     generators - but with the extra requirement of having
//                     a non-zero electric charge.
//                     While for G4UnknownParticle is possible to assign
//                     transportation and decay processes, for
//                     G4ChargedUnknownParticle is also possible to assign
//                     ionisation and multiple scattering.
//
//      Modifications:
//      
// --------------------------------------------------------------------------
//

#ifndef G4ChargedUnknownParticle_h
#define G4ChargedUnknownParticle_h 1

#include "G4ParticleDefinition.hh"


class G4ChargedUnknownParticle : public G4ParticleDefinition {
  public:
    static G4ChargedUnknownParticle* Definition();
    static G4ChargedUnknownParticle* ChargedUnknownParticleDefinition();
    static G4ChargedUnknownParticle* ChargedUnknownParticle();
  private:
  G4ChargedUnknownParticle() = default;
    ~G4ChargedUnknownParticle() = default;
    static G4ChargedUnknownParticle* theInstance;
};

#endif
