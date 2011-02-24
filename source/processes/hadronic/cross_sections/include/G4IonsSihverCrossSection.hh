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
#ifndef G4IonsSihverCrossSection_h
#define G4IonsSihverCrossSection_h
//
// Class Description
// Implementation of formulas 
// Sihver et al. Phys. Rev. C 47 1225 (1993); 
// Total Reaction Cross Section for Nucleus-nucles reactions.
//    Energy independent   
//    Valid for 100>MeV/nucleon 
// Class Description - End
//
// 23-Dec-2006 Isotope dependence added by D. Wright
//

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4IonsSihverCrossSection : public G4VCrossSectionDataSet
{
  public:
    G4IonsSihverCrossSection()
     : G4VCrossSectionDataSet("G4IonsSihverCrossSection"),
       square_r0 ( (1.36*fermi) * (1.36*fermi) )
    {}
   
    virtual
    G4bool IsApplicable(const G4DynamicParticle* aDP, const G4Element*)
    {
      return IsIsoApplicable(aDP, 0, 0);
    }

    virtual 
    G4bool IsIsoApplicable(const G4DynamicParticle* aDP, G4int /*ZZ*/, 
                           G4int /*AA*/)
    {
      G4int BaryonNumber = aDP->GetDefinition()->GetBaryonNumber();
      G4double KineticEnergy = aDP->GetKineticEnergy(); 
      if ( KineticEnergy / BaryonNumber >= 100*MeV && BaryonNumber > 1 ) 
         return true;
      return false;
    }

    virtual
    G4double GetCrossSection(const G4DynamicParticle*, 
                             const G4Element*, G4double aTemperature);

    virtual
    G4double GetZandACrossSection(const G4DynamicParticle*, G4int ZZ, 
                                  G4int AA, G4double aTemperature);

    virtual
    void BuildPhysicsTable(const G4ParticleDefinition&)
    {}

    virtual
    void DumpPhysicsTable(const G4ParticleDefinition&) 
    {G4cout << "tG4GIonCrossSection: uses formula"<<G4endl;}

  private:
    const G4double square_r0;

};

#endif
