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
// 31 Jul 2001   MGP                 Created
//
// 15 Jul 2009   N.A.Karakatsanis    New methods added for loading logarithmic data
//                                   to enhance computing performance of interpolation
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// -------------------------------------------------------------------

#ifndef G4EMDATASET_HH
#define G4EMDATASET_HH 1

#include <CLHEP/Units/SystemOfUnits.h>
#include "globals.hh"
#include "G4VEMDataSet.hh"

class G4VDataSetAlgorithm;

class G4EMDataSet : public G4VEMDataSet
{
public:
  explicit G4EMDataSet(G4int argZ, 
		       G4VDataSetAlgorithm* algo, 
		       G4double xUnit=CLHEP::MeV, 
		       G4double yUnit=CLHEP::barn,
		       G4bool random=false);

  explicit G4EMDataSet(G4int argZ, 
		       G4DataVector* xData, 
		       G4DataVector* data, 
		       G4VDataSetAlgorithm* algo, 
		       G4double xUnit=CLHEP::MeV, 
		       G4double yUnit=CLHEP::barn,
		       G4bool random=false);

  explicit G4EMDataSet(G4int argZ, 
		       G4DataVector* xData, 
		       G4DataVector* data,
		       G4DataVector* xLogData, 
		       G4DataVector* Logdata, 
		       G4VDataSetAlgorithm* algo, 
		       G4double xUnit=CLHEP::MeV, 
		       G4double yUnit=CLHEP::barn,
		       G4bool random=false);

  virtual ~G4EMDataSet();
 
  virtual G4double FindValue(G4double x, G4int componentId=0) const;
  
  virtual void PrintData(void) const;

  virtual const G4VEMDataSet* GetComponent(G4int /* componentId */) const { return 0; }

  virtual void AddComponent(G4VEMDataSet* /* dataSet */) {}

  virtual std::size_t NumberOfComponents(void) const { return 0; }

  virtual const G4DataVector& GetEnergies(G4int /* componentId */) const { return *energies; }
  virtual const G4DataVector& GetData(G4int /* componentId */) const { return *data; }
  virtual const G4DataVector& GetLogEnergies(G4int /* componentId */) const { return *log_energies; }
  virtual const G4DataVector& GetLogData(G4int /* componentId */) const { return *log_data; }

  virtual void SetEnergiesData(G4DataVector* xData, G4DataVector* data, G4int componentId);
  virtual void SetLogEnergiesData(G4DataVector* xData,
                                  G4DataVector* data,
                                  G4DataVector* xLogData, 
                                  G4DataVector* Logdata,
                                  G4int componentId);


  virtual G4bool LoadData(const G4String& fileName);
  virtual G4bool LoadNonLogData(const G4String& fileName);

  virtual G4bool SaveData(const G4String& fileName) const;

  virtual G4double RandomSelect(G4int componentId = 0) const;
    

private:

  std::size_t FindLowerBound(G4double energy) const;
  std::size_t FindLowerBound(G4double x, G4DataVector* values) const;

  G4double IntegrationFunction(G4double x);
  virtual void BuildPdf();
  G4String FullFileName(const G4String& fileName) const;

  // Hide copy constructor and assignment operator 
  explicit G4EMDataSet();
  G4EMDataSet(const G4EMDataSet& copy) = delete;
  G4EMDataSet& operator=(const G4EMDataSet& right) = delete;

  G4DataVector* energies;            // Owned pointer
  G4DataVector* data;                // Owned pointer
  G4DataVector* log_energies;        // Owned pointer
  G4DataVector* log_data;            // Owned pointer
  G4VDataSetAlgorithm* algorithm;    // Owned pointer 
  G4DataVector* pdf;

  G4double unitEnergies;
  G4double unitData;

  G4int z;
  G4bool randomSet;
};
#endif /* G4EMDATASET_HH */
