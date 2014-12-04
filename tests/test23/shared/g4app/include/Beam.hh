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
#ifndef Beam_H
#define Beam_H 1

#include "G4LorentzVector.hh"
#include "G4String.hh"

class Beam
{

   public:
   
      Beam() : fBeamPartName(""), fBeamPartMass(0.), fBeamEnergy(0.),
                fLabV(G4LorentzVector()), fLabP(G4LorentzVector()) {}
      ~Beam() {}
      
      void SetBeam( G4String name, G4double mass, G4double energy ) { fBeamPartName=name;
                                                                      fBeamPartMass=mass;
						                      fBeamEnergy=energy;
						                      return; }
						 
      void SetLabV( G4LorentzVector  lv ) { fLabV=lv; return; }
      void SetLabP( G4LorentzVector  lp ) { fLabP=lp; return; }
      
      G4String         GetBeamPartName()   const { return fBeamPartName; }
      G4double         GetBeamPartMass()   const { return fBeamPartMass; }
      G4double         GetBeamEnergy()     const { return fBeamEnergy; }
      const G4LorentzVector&  GetLabV()    const { return fLabV; }
      const G4LorentzVector&  GetLabP()    const { return fLabP; }

   private:
   
      G4String        fBeamPartName;
      G4double        fBeamPartMass;
      G4double        fBeamEnergy;
      G4LorentzVector fLabV;
      G4LorentzVector fLabP;

};

#endif
