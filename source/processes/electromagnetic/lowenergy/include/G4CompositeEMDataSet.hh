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
// $Id: G4CompositeEMDataSet.hh,v 1.4 2002-05-28 09:15:26 pia Exp $
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
// Composite data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4COMPOSITEEMDATASET_HH
#define G4COMPOSITEEMDATASET_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4VEMDataSet.hh"
#include "g4std/vector"

class G4EMDataSet;
class G4VDataSetAlgorithm;

class G4CompositeEMDataSet : public G4VEMDataSet 
{ 
public:

  G4CompositeEMDataSet(G4VDataSetAlgorithm* interpolation,
		       G4double unitE = MeV, G4double unitData = barn,
		       G4int minZ = 1, G4int maxZ = 99);

  G4CompositeEMDataSet(const G4String& dataFile,
		       G4VDataSetAlgorithm* interpolation,
		       G4double unitE = MeV, G4double unitData = barn,
		       G4int minZ = 1, G4int maxZ = 99);

  ~G4CompositeEMDataSet();
 
  G4double FindValue(G4double e, G4int Z = 0) const;

  const G4VEMDataSet* GetComponent(G4int i) const { return components[i]; }

  void AddComponent(G4VEMDataSet* component);
  
  size_t NumberOfComponents() const { return nComponents; }

  void PrintData() const;

  const G4DataVector& GetEnergies(G4int i) const;
  const G4DataVector& GetData(G4int i) const;

private:

  // Hide copy constructor and assignment operator 
  G4CompositeEMDataSet& operator=(const G4CompositeEMDataSet& right);
  G4CompositeEMDataSet(const G4CompositeEMDataSet&);

  void LoadData(const G4String& fileName);

  G4VDataSetAlgorithm* algorithm; 

  G4std::vector<G4VEMDataSet*> components;
  size_t nComponents;

  G4double unit1;
  G4double unit2;

  G4int zMin;
  G4int zMax;
};
 
#endif









