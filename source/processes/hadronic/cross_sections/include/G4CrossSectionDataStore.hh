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
                            const G4Element*, G4double aTemperature);

   void AddDataSet(G4VCrossSectionDataSet*);

   void BuildPhysicsTable(const G4ParticleDefinition&);

   void DumpPhysicsTable(const G4ParticleDefinition&);

   void SetVerboseLevel(G4int value)
   {
      verboseLevel = value;
   }

   G4int GetVerboseLevel()
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
