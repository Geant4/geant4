// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronInelasticDataSet.hh,v 1.2 1999-12-15 14:52:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4HadronInelasticDataSet -- header file
// F.W. Jones, TRIUMF, 19-MAY-98
//

#ifndef G4HadronInelasticDataSet_h
#define G4HadronInelasticDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4HadronInelasticDataSet : public G4VCrossSectionDataSet
{
public:

   G4HadronInelasticDataSet()
   {
      theHadronCrossSections = G4HadronCrossSections::Instance();
   }

   ~G4HadronInelasticDataSet()
   {
   }

   G4bool IsApplicable(const G4DynamicParticle* aParticle,
                       const G4Element* anElement)
   {
      return theHadronCrossSections->IsApplicable(aParticle, anElement);
   }

   G4double GetCrossSection(const G4DynamicParticle* aParticle,
                            const G4Element* anElement)
   {
      return theHadronCrossSections->GetInelasticCrossSection(aParticle,
                                                              anElement);
   }

   void BuildPhysicsTable(const G4ParticleDefinition&)
   {
   }

   void DumpPhysicsTable(const G4ParticleDefinition&)
   {
   }

private:

   G4HadronCrossSections* theHadronCrossSections;
};

#endif
