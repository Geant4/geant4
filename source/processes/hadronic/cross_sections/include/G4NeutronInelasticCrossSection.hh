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

#ifndef G4NeutronInelasticCrossSection_h
#define G4NeutronInelasticCrossSection_h

#include "globals.hh"
#include "G4Neutron.hh"

// Class Description
// Cross-sections for neutron nuclear scattering from 14 MeV up to 20 GeV, getting the
// low energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), M. Laidlaw (British Columbia U.). 1996. 
// Class Description - End

#include "G4VCrossSectionDataSet.hh"

class G4NeutronInelasticCrossSection : public G4VCrossSectionDataSet
{
   public:
   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element* aEle)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()==G4Neutron::Neutron()) &&
        ( aPart->GetKineticEnergy()<20*GeV) &&
          aPart->GetKineticEnergy()>19.9*MeV) result = true;
     if(aEle->GetZ()<2) result = false;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, 
                            const G4Element*, G4double aTemperature);
   
   G4double GetCrossSection(G4double anEnergy, G4double anA, G4double aZ);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4NeutronInelasticCrossSection: uses formula"<<G4endl;}

};

#endif
