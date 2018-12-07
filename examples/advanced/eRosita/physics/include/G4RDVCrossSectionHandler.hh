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
// 16 Sep 2001   MGP           Created
// 26 Sep 2001   V.Ivanchenko  Hide copy constructor and assignement operator
// 18 Apr 2002   V.Ivanchenko  Move member function ValueForMaterial to public
// 21 Jan 2003   V.Ivanchenko  Cut per region
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Base class for cross section manager for an electromagnetic physics process
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDVCROSSSECTIONHANDLER_HH
#define G4RDVCROSSSECTIONHANDLER_HH 1

#include <map>
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4MaterialCutsCouple.hh"

class G4RDVDataSetAlgorithm;
class G4RDVEMDataSet;
class G4Material;
class G4Element;

class G4RDVCrossSectionHandler {

public:

  G4RDVCrossSectionHandler();

  G4RDVCrossSectionHandler(G4RDVDataSetAlgorithm* interpolation,
			 G4double minE = 250*CLHEP::eV, G4double maxE = 100*CLHEP::GeV,
			 G4int nBins = 200,
			 G4double unitE = CLHEP::MeV, G4double unitData = CLHEP::barn,
			 G4int minZ = 1, G4int maxZ = 99);

  virtual ~G4RDVCrossSectionHandler();

  void Initialise(G4RDVDataSetAlgorithm* interpolation = 0,
		  G4double minE = 250*CLHEP::eV, G4double maxE = 100*CLHEP::GeV,
		  G4int numberOfBins = 200,
		  G4double unitE = CLHEP::MeV, G4double unitData = CLHEP::barn,
		  G4int minZ = 1, G4int maxZ = 99);

  G4int SelectRandomAtom(const G4MaterialCutsCouple* couple, G4double e) const;

  const G4Element* SelectRandomElement(const G4MaterialCutsCouple* material,
				             G4double e) const;

  G4int SelectRandomShell(G4int Z, G4double e) const;

  G4RDVEMDataSet* BuildMeanFreePathForMaterials(const G4DataVector* energyCuts = 0);

  G4double FindValue(G4int Z, G4double e) const;

  G4double FindValue(G4int Z, G4double e, G4int shellIndex) const;

  G4double ValueForMaterial(const G4Material* material, G4double e) const;

  void LoadData(const G4String& dataFile);

  void LoadShellData(const G4String& dataFile);

  void PrintData() const;

  void Clear();

protected:

  G4int NumberOfComponents(G4int Z) const;

  void ActiveElements();

  // Factory method
  virtual std::vector<G4RDVEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
								       const G4DataVector* energyCuts = 0) = 0;

  // Factory method
  virtual G4RDVDataSetAlgorithm* CreateInterpolation();

  const G4RDVDataSetAlgorithm* GetInterpolation() const { return interpolation; }


private:

  // Hide copy constructor and assignment operator
  G4RDVCrossSectionHandler(const G4RDVCrossSectionHandler&);
  G4RDVCrossSectionHandler & operator=(const G4RDVCrossSectionHandler &right);

  G4RDVDataSetAlgorithm* interpolation;

  G4double eMin;
  G4double eMax;
  G4int nBins;

  G4double unit1;
  G4double unit2;

  G4int zMin;
  G4int zMax;

  G4DataVector activeZ;

  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> > dataMap;

  std::vector<G4RDVEMDataSet*>* crossSections;
};

#endif
