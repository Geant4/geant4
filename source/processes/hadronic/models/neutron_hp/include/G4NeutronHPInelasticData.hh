// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPInelasticData.hh,v 1.3 2000-12-14 09:20:35 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPInelasticData_h
#define G4NeutronHPInelasticData_h 1

// Class Description
// Cross-section data set for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"

class G4NeutronHPInelasticData : public G4VCrossSectionDataSet
{
   public:
   
   G4NeutronHPInelasticData();
   
   ~G4NeutronHPInelasticData();
   
   G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*);

   void BuildPhysicsTable(const G4ParticleDefinition&);

   void DumpPhysicsTable(const G4ParticleDefinition&);
   
   private:
   
   G4PhysicsTable * theCrossSections;
};

#endif
