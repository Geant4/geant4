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
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDEMDATASET_HH
#define G4RDEMDATASET_HH 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4RDVEMDataSet.hh"

class G4RDVDataSetAlgorithm;

class G4RDEMDataSet : public G4RDVEMDataSet
{
public:
  G4RDEMDataSet(G4int argZ, 
	      G4RDVDataSetAlgorithm* algo, 
	      G4double xUnit=CLHEP::MeV, 
	      G4double yUnit=CLHEP::barn,
              G4bool random=false);

  G4RDEMDataSet(G4int argZ, 
	      G4DataVector* xData, 
	      G4DataVector* data, 
	      G4RDVDataSetAlgorithm* algo, 
	      G4double xUnit=CLHEP::MeV, 
	      G4double yUnit=CLHEP::barn,
	      G4bool random=false);

  virtual ~G4RDEMDataSet();
 
  virtual G4double FindValue(G4double x, G4int componentId=0) const;
  
  virtual void PrintData(void) const;

  virtual const G4RDVEMDataSet* GetComponent(G4int /* componentId */) const { return 0; }

  virtual void AddComponent(G4RDVEMDataSet* /* dataSet */) {}

  virtual size_t NumberOfComponents(void) const { return 0; }

  virtual const G4DataVector& GetEnergies(G4int /* componentId */) const { return *energies; }
  virtual const G4DataVector& GetData(G4int /* componentId */) const { return *data; }
  virtual void SetEnergiesData(G4DataVector* xData, G4DataVector* data, G4int componentId);

  virtual G4bool LoadData(const G4String& fileName);
  virtual G4bool SaveData(const G4String& fileName) const;

  virtual G4double RandomSelect(G4int componentId = 0) const;
    

private:

  size_t FindLowerBound(G4double energy) const;
  size_t FindLowerBound(G4double x, G4DataVector* values) const;

  G4double IntegrationFunction(G4double x);

  virtual void BuildPdf();
  
  G4String FullFileName(const G4String& fileName) const;

  // Hide copy constructor and assignment operator 
  G4RDEMDataSet();
  G4RDEMDataSet(const G4RDEMDataSet& copy);
  G4RDEMDataSet& operator=(const G4RDEMDataSet& right);

  G4int z;

  G4DataVector* energies;            // Owned pointer
  G4DataVector* data;                // Owned pointer

  G4RDVDataSetAlgorithm* algorithm;    // Owned pointer 
  
  G4double unitEnergies;
  G4double unitData;

  G4DataVector* pdf;
  G4bool randomSet;
};
#endif /* G4EMDATASET_HH */
