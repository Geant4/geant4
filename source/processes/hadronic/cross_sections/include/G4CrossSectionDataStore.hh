// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CrossSectionDataStore.hh,v 1.3 2000-12-14 08:33:02 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4CrossSectionDataStore -- header file
// F.W. Jones, TRIUMF, 19-NOV-97
//
// Class Description
// This is the class to which to register data-sets. You can get the instance
// from energy hadronic process, and use its 'AddDataSet(...)' method to tailor
// the cross-sectinos for your application.
// Class Description - End

#ifndef G4CrossSectionDataStore_h
#define G4CrossSectionDataStore_h 1

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4VCrossSectionDataSet.hh"


class G4CrossSectionDataStore
{
public:

   G4CrossSectionDataStore() :
      NDataSetList(0), verboseLevel(0)
   {
   }

   ~G4CrossSectionDataStore()
   {
   }

   G4double GetCrossSection(const G4DynamicParticle*, 
                            const G4Element*);

   void AddDataSet(G4VCrossSectionDataSet*);

   void BuildPhysicsTable(const G4ParticleDefinition&);

   void DumpPhysicsTable(const G4ParticleDefinition&);

   void SetVerboseLevel(G4int value)
   {
      verboseLevel = value;
   }

   G4int GetVerboseLevel(G4int value)
   {
      return verboseLevel;
   }

private:

   enum { NDataSetMax = 100 };
   G4VCrossSectionDataSet* DataSetList[NDataSetMax];
   G4int NDataSetList;
   G4int verboseLevel;
};

#endif
