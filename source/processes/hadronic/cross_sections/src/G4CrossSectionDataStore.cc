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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// GEANT4 physics class: G4CrossSectionDataStore
// F.W. Jones, TRIUMF, 19-NOV-97
//

#include "G4CrossSectionDataStore.hh"
#include "G4HadronicException.hh"


G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         const G4Element* anElement,
					 G4double aTemperature)
{
   if (NDataSetList == 0) 
   {
      throw G4HadronicException(__FILE__, __LINE__, 
       "G4CrossSectionDataStore: no data sets registered");
      return DBL_MIN;
   }
   for (G4int i = NDataSetList-1; i >= 0; i--) {
      if (DataSetList[i]->IsApplicable(aParticle, anElement))
             return DataSetList[i]->GetCrossSection(aParticle, anElement, aTemperature);
   }
   throw G4HadronicException(__FILE__, __LINE__, 
                            "G4CrossSectionDataStore: no applicable data set found "
                            "for particle/element");
   return DBL_MIN;
}


void
G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* aDataSet)
{
   if (NDataSetList == NDataSetMax) {
      G4cout << "WARNING: G4CrossSectionDataStore::AddDataSet: "<<G4endl;
      G4cout << "         reached maximum number of data sets";
      G4cout << "         data set not added !!!!!!!!!!!!!!!!";
      return;
   }
   DataSetList[NDataSetList] = aDataSet;
   NDataSetList++;
}


void
G4CrossSectionDataStore::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (NDataSetList == 0) 
   {
     G4Exception("G4CrossSectionDataStore", "007", FatalException,
                 "BuildPhysicsTable: no data sets registered");
     return;
   }
   for (G4int i = NDataSetList-1; i >= 0; i--) {
      DataSetList[i]->BuildPhysicsTable(aParticleType);
   }
}


void
G4CrossSectionDataStore::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
   if (NDataSetList == 0) {
      G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: no data sets registered"<<G4endl;
      return;
   }
   for (G4int i = NDataSetList-1; i >= 0; i--) {
      DataSetList[i]->DumpPhysicsTable(aParticleType);
   }
}
