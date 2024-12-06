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
//---------------------------------------------------------------------------
//
// ClassName:    G4URRNeutrons
//
// Author:       Alberto Ribon - October 2024  
//
// Description:  Physics list constructor that can be applied on top of any
//               _HP or _HPT based physics list.
//               This class enables the special Unresolved Resonance Region
//               (URR) treatment of low-energy neutrons based on Particle
//               Table (PT).
//               If this constructor is applied on top of a non-HP based
//               physics list, then nothing is done (i.e. the physics list
//               remains as it was originally, and a warning is printed out).
//
// Modified:
//
//----------------------------------------------------------------------------
//
// Addition of neutron thermal scattering on top of any PhysicsList

#ifndef G4URRNeutrons_h
#define G4URRNeutrons_h

#include "G4VHadronPhysics.hh"
#include "globals.hh"


class G4URRNeutrons : public G4VHadronPhysics {
  public:
    explicit G4URRNeutrons( G4int ver = 1 );
    virtual ~G4URRNeutrons();
    void ConstructProcess() override;
};

#endif
