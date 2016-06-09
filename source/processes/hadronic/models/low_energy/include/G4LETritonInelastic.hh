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
// $Id$
//
// Hadronic Process: Low Energy Triton Inelastic Process
// J.L. Chuma, TRIUMF, 25-Feb-1997
 
#ifndef G4LETritonInelastic_h
#define G4LETritonInelastic_h 1
 
// Class Description
// Final state production model for Triton inelastic scattering below 100 MeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4InelasticInteraction.hh"
 
class G4LETritonInelastic : public G4InelasticInteraction
{
  public:
    
    G4LETritonInelastic() : G4InelasticInteraction("G4LETritonInelastic")
    {
      SetMinEnergy( 0.0 );
      // SetMaxEnergy( 100.*CLHEP::MeV );  // NUCREC only worked for energies < 100MeV
      // Work around to avoid exception in G4EnergyRangeManager
      SetMaxEnergy(10.*CLHEP::TeV);  // NUCREC only worked for energies < 100MeV
      G4cout << "WARNING: model G4LETritonInelastic is being deprecated and will\n"
             << "disappear in Geant4 version 10.0"  << G4endl;
    }
    
    ~G4LETritonInelastic() { }
    
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& targetNucleus);

    virtual void ModelDescription(std::ostream& outFile) const;
 };
 
#endif
 
