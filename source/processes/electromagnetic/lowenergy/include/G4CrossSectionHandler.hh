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
// $Id: G4CrossSectionHandler.hh,v 1.3 2001-09-05 12:29:49 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
//  1 Aug 2001   MGP        Created
//
//  Modified: 30.08.01 V.Ivanchenko add G4VEMSecondaryGenerator
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set manager for an electromagnetic physics process
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONHANDLER_HH
#define G4CROSSSECTIONHANDLER_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "g4std/map"
#include "g4std/vector"

class G4VDataSetAlgorithm;
class G4VEMDataSet;
class G4Material;
class G4Element;
class G4VEMSecondaryGenerator;

class G4CrossSectionHandler {
 
public:

  G4CrossSectionHandler(const G4VDataSetAlgorithm* interpolation,
			G4double minE = 250*eV, G4double maxE = 100*GeV, 
			G4int nBins = 200,
			G4double unitE = MeV, G4double unitData = barn,
			G4int minZ = 1, G4int maxZ = 99);

  ~G4CrossSectionHandler();
	
  G4double FindValue(G4int Z, G4double e) const;

  G4double ValueForMaterial(const G4Material* material, G4double e) const;

  G4int SelectRandomAtom(const G4Material* material, G4double e) const;

  const G4Element* SelectRandomElement(const G4Material* material, G4double e) const;

  G4int SelectRandomShell(G4int Z, G4double e) const;

  // For Discrete Processes
  G4VEMDataSet* BuildMeanFreePathForMaterials();

  G4VEMDataSet* BuildMeanFreePathForMaterials(
                const G4DataVector& minEnergy,
                const G4DataVector& maxEnergy);
 
  void LoadData(const G4String& dataFile);

  void LoadShellData(const G4String& dataFile);

  void PrintData() const;
  
  void Clear();

  void SetSecondaryGenerator(G4VEMSecondaryGenerator* p) {theGenerator = p;};

  void SetVerbose(G4int val) {verbose = val;};
   
private:
 
  // Hide copy constructor and assignment operator

  void BuildCrossSectionsForMaterials(const G4DataVector& energyVector);

  void BuildCrossSectionsWithCut(const G4DataVector& energyVector,
                                 const G4DataVector& minEnergy,
                                 const G4DataVector& maxEnergy);

  G4VEMDataSet* BuildMeanFreePathForMaterials(const G4DataVector& energies);

  void ActiveElements();

  const G4VDataSetAlgorithm* interpolation;
  
  G4double eMin;
  G4double eMax;
  G4int nBins;

  G4double unit1;
  G4double unit2;

  G4int zMin;
  G4int zMax; 

  G4DataVector activeZ;

  // Data from database
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> > dataMap;

  // Calculated data
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> > dataWithCutMap;
  
  G4std::vector<G4VEMDataSet*> crossSections;

  G4VEMSecondaryGenerator* theGenerator;
 
  G4int verbose;
};
 
#endif











