// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronCaptureDataSet.hh,v 1.1 1999-01-07 16:11:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4HadronCaptureDataSet -- header file
// F.W. Jones, TRIUMF, 19-MAY-98
//

#ifndef G4HadronCaptureDataSet_h
#define G4HadronCaptureDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4HadronCaptureDataSet : public G4VCrossSectionDataSet
{
public:

   G4HadronCaptureDataSet()
   {
      theHadronCrossSections = G4HadronCrossSections::Instance();
   }

   ~G4HadronCaptureDataSet()
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
      return theHadronCrossSections->GetCaptureCrossSection(aParticle,
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
