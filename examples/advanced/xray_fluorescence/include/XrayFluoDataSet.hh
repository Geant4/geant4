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
// $Id: XrayFluoDataSet.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------

//for the documentation of this class see also G4EMDataSet class

#ifndef FluoDataSet_hh
#define FluoDataSet_hh 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4VEMDataSet.hh"

class G4VDataSetAlgorithm;

class XrayFluoDataSet : public G4VEMDataSet {

public:

  XrayFluoDataSet(G4int Z,
		  G4DataVector* points, 
	      G4DataVector* values,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);
  
  XrayFluoDataSet(G4int Z,
		  const G4String& dataFile,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  ~XrayFluoDataSet();
 
  //find the value corresponding to the energy e in the set
  //identified by id
  G4double FindValue(G4double e, G4int id = 0) const;

  virtual const G4VEMDataSet* GetComponent(G4int i) const { return 0;} 

  virtual void AddComponent(G4VEMDataSet* dataSet) { }

  virtual size_t NumberOfComponents() const { return 0; }


  void PrintData() const;

  const G4DataVector& GetEnergies(G4int i) const { return *energies; }
  const G4DataVector& GetData(G4int i) const { return *data; }

private:

  void LoadData(const G4String& dataFile);
  G4int z;
  G4int FindBinLocation(G4double energy) const;

  G4DataVector* energies; // Owned pointer
  G4DataVector* data;     // Owned pointer

  const G4VDataSetAlgorithm* algorithm; // Not owned pointer 
  
  G4double unit1;
  G4double unit2;

  size_t numberOfBins;

};
 
#endif
