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
// by JPW, working, but to be cleaned up. @@@@

#ifndef G4TripathiCrossSection_h
#define G4TripathiCrossSection_h
// Class Description
// Implementation of formulas in analogy to NASA technical paper 3621 by 
// Tripathi, et al.; Cross-sections for ion ion scattering.

// Class Description - End

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4TripathiCrossSection : public G4VCrossSectionDataSet
{
   public:
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element*)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()->GetBaryonNumber()>2.5) &&
        ( aPart->GetKineticEnergy()/aPart->GetDefinition()->GetBaryonNumber()<1*GeV) ) result = true;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, 
                            const G4Element*, G4double aTemperature);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4TripathiCrossSection: uses formula"<<G4endl;}

};

#endif
