// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronCrossSectionPlugin.hh,v 1.1 1999-01-08 16:33:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Plug-in for G4CrossSectionDataTest
// F.W. Jones, TRIUMF, 30-MAR-98
//  


#ifndef G4HadronCrossSectionPlugin_h
#define G4HadronCrossSectionPlugin_h 1
 
#include "G4VCrossSectionDataSet.hh"

#include "globals.hh"
#include "G4Element.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"


class G4HadronCrossSectionPlugin : public G4VCrossSectionDataSet
{
public:

   G4HadronCrossSectionPlugin()
   {
   }

   ~G4HadronCrossSectionPlugin()
   {
   }

   G4double
   GetCrossSection(const G4DynamicParticle* aParticle,
                   const G4Element* anElement)
   {
      return 999.*millibarn;
   }

   G4bool
   IsApplicable(const G4DynamicParticle* aParticle,
                const G4Element* anElement)
   {
      if (verboseLevel > 1) {
         G4cout << "G4HadronCrossSectionPlugin::IsApplicable:" << endl;
         G4cout << "  Particle: " <<
                 aParticle->GetDefinition()->GetParticleName() << endl;
         G4cout << "  Energy:   " << aParticle->GetKineticEnergy()/GeV << endl;
         G4cout << "  Element:  " << anElement->GetName() << endl;
      }
      if (aParticle->GetDefinition() != G4Proton::Proton()) return 0;
      G4double ekin = aParticle->GetKineticEnergy()/GeV;
      return (ekin > 0.95 && ekin < 1.05);
   }

   void BuildPhysicsTable(const G4ParticleDefinition&)
   {
   }

   void DumpPhysicsTable(const G4ParticleDefinition&)
   {
   }

private:

};
#endif
