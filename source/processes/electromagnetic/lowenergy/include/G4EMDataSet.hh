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
// $Id: G4EMDataSet.hh,v 1.1 2001-08-20 16:36:01 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef G4EMDATASET_HH
#define G4EMDATASET_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4VEMDataSet.hh"

class G4VDataSetAlgorithm;

class G4EMDataSet : public G4VEMDataSet {
 
public:

  G4EMDataSet(G4int Z, 
	      G4DataVector* points, 
	      G4DataVector* values,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  G4EMDataSet(G4int Z, 
	      const G4String& dataFile,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  ~G4EMDataSet();
 
  G4double FindValue(G4double e, G4int id = 0) const;
  
  void PrintData() const;

  const G4DataVector& GetEnergies(G4int i) const { return *energies; }
  const G4DataVector& GetData(G4int i) const { return *data; }

private:

  void LoadData(const G4String& dataFile);
  
  G4int FindBinLocation(G4double energy) const;

  G4int z;

  G4DataVector* energies; // Owned pointer
  G4DataVector* data;     // Owned pointer

  const G4VDataSetAlgorithm* algorithm; // Not owned pointer 
  
  G4double unit1;
  G4double unit2;

  size_t numberOfBins;

};
 
#endif
 










