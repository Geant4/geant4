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
#ifndef G4IonProtonCrossSection_h
#define G4IonProtonCrossSection_h

#include "globals.hh"
#include "G4Proton.hh"
// Class Description
// Cross-sections for ion proton scattering up to 20 GeV, getting the low
// energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), D. Axen (British Columbia U.). 1996. 
// Published in Phys.Rev.C54:1329-1332,1996 
// Class Description - End

#include "G4VCrossSectionDataSet.hh"
#include "G4ProtonInelasticCrossSection.hh"

class G4IonProtonCrossSection : public G4VCrossSectionDataSet
{
   public:

   virtual
   G4bool IsZAApplicable(const G4DynamicParticle* aPart, G4double /*ZZ*/,
                         G4double AA)
   {
     G4bool result = false;
     if((AA < 1.1) &&
        ( aPart->GetKineticEnergy()/aPart->GetDefinition()
                                         ->GetBaryonNumber() < 20*GeV &&
	  aPart->GetDefinition()->GetBaryonNumber() > 4)
       ) result = true;
     return result;
   }

   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element* anEle)
   {
     return IsZAApplicable(aPart, 0., anEle->GetN());
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle* aPart, 
                            const G4Element*, G4double )
   {
     return GetIsoZACrossSection(aPart, 0., 0., 0.);
   }

   virtual
   G4double GetIsoZACrossSection(const G4DynamicParticle* aPart, 
                                 G4double /*ZZ*/, G4double /*AA*/, 
                                 G4double /*temperature*/)
   {
     G4ProtonInelasticCrossSection theForward;
     return theForward.GetCrossSection(aPart->GetKineticEnergy(),
                                  aPart->GetDefinition()->GetBaryonNumber(),
				  aPart->GetDefinition()->GetPDGCharge());
   }


   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&)
   {G4cout << "G4IonProtonCrossSection: uses formula"<<G4endl;}

};

#endif
