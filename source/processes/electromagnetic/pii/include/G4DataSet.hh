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
// $Id: G4DataSet.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created as G4EMDataSet
// 31 Jul 2008   MGP        Revised and renamed to G4DataSet
//
// -------------------------------------------------------------------

// Class description:
// Leaf data set 
// Applies a Composite design pattern for data library management

// -------------------------------------------------------------------

#ifndef  G4DATASET_HH
#define  G4DATASET_HH 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4IDataSet.hh"

class G4IInterpolator;

class G4DataSet : public G4IDataSet
{
public:
  G4DataSet(G4int argZ, 
	    G4IInterpolator* algo, 
	    G4double xUnit=CLHEP::MeV, 
	    G4double yUnit=CLHEP::barn,
	    G4bool random=false);

  G4DataSet(G4int argZ, 
	    G4DataVector* xData, 
	    G4DataVector* data, 
	    G4IInterpolator* algo, 
	    G4double xUnit=CLHEP::MeV, 
	    G4double yUnit=CLHEP::barn,
	    G4bool random=false);

  virtual ~G4DataSet();
 
  virtual G4double FindValue(G4double x, G4int componentId=0) const;
  
  virtual void PrintData(void) const;

  virtual const G4IDataSet* GetComponent(G4int /* componentId */) const { return 0; }

  virtual void AddComponent(G4IDataSet* /* dataSet */) {}

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
  G4DataSet();
  G4DataSet(const G4DataSet& copy);
  G4DataSet& operator=(const G4DataSet& right);

  G4int z;

  G4DataVector* energies;            // Owned pointer
  G4DataVector* data;                // Owned pointer

  G4IInterpolator* algorithm;    // Owned pointer 
  
  G4double unitEnergies;
  G4double unitData;

  G4DataVector* pdf;
  G4bool randomSet;
};
#endif /* G4DATASET_HH */
