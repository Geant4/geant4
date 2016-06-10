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
// $Id: G4PixeCrossSectionHandler.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 16 Jun 2008   MGP           Created on the basis of G4CrossSectionHandler

// -------------------------------------------------------------------
// Class description:
// Cross section manager for hadron impact ionization
// Documented in:
// M.G. Pia et al., PIXE Simulation With Geant4,
// IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, Dec. 2009.

// -------------------------------------------------------------------

#ifndef G4PIXECROSSSECTIONHANDLER_HH
#define G4PIXECROSSSECTIONHANDLER_HH 1

#include <map>
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4MaterialCutsCouple.hh"

class G4IInterpolator;
class G4IDataSet;
class G4Material;
class G4Element;

class G4PixeCrossSectionHandler {

public:

  G4PixeCrossSectionHandler();

  G4PixeCrossSectionHandler(G4IInterpolator* interpolation,
			    const G4String& modelK="ecpssr",
			    const G4String& modelL="ecpssr",
			    const G4String& modelM="ecpssr",
			    G4double minE = 1*CLHEP::keV,
                            G4double maxE = 0.1*CLHEP::GeV,
			    G4int nBins = 200,
			    G4double unitE = CLHEP::MeV,
                            G4double unitData = CLHEP::barn,
			    G4int minZ = 6, G4int maxZ = 92);
  
  virtual ~G4PixeCrossSectionHandler();

  void Initialise(G4IInterpolator* interpolation,
		  const G4String& modelK="ecpssr",
		  const G4String& modelL="ecpssr",
		  const G4String& modelM="ecpssr",
		  G4double minE = 1*CLHEP::keV,
                  G4double maxE = 0.1*CLHEP::GeV,
		  G4int nBins = 200,
		  G4double unitE = CLHEP::MeV,
                  G4double unitData = CLHEP::barn,
		  G4int minZ = 6, G4int maxZ = 92);

  G4int SelectRandomAtom(const G4Material* material, G4double e) const;

  G4int SelectRandomShell(G4int Z, G4double e) const;

  G4double FindValue(G4int Z, G4double e) const;

  G4double FindValue(G4int Z, G4double e, G4int shellIndex) const;

  G4double ValueForMaterial(const G4Material* material, G4double e) const;

  // void LoadData(const G4String& dataFile);

  void LoadShellData(const G4String& dataFile);

  // Ionisation cross section as in Geant4 Physics Reference Manual
  G4double MicroscopicCrossSection(const G4ParticleDefinition* particleDef,
				   G4double kineticEnergy,
				   G4double Z,
				   G4double deltaCut) const;

  void PrintData() const;

  void Clear();

private:

 // Hide copy constructor and assignment operator
  G4PixeCrossSectionHandler(const G4PixeCrossSectionHandler&);
  G4PixeCrossSectionHandler & operator=(const G4PixeCrossSectionHandler &right);

  G4int NumberOfComponents(G4int Z) const;

  void ActiveElements();

  void BuildForMaterials();
  // Factory method
  std::vector<G4IDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector);

  // Factory method
  G4IInterpolator* CreateInterpolation();

  const G4IInterpolator* GetInterpolation() const { return interpolation; }
 
  G4IInterpolator* interpolation;

  G4double eMin;
  G4double eMax;
  G4int nBins;

  G4double unit1;
  G4double unit2;

  G4int zMin;
  G4int zMax;

  G4DataVector activeZ;

  // Map of PixeShellDataSets with the shell cross sections for each element
  std::map<G4int,G4IDataSet*,std::less<G4int> > dataMap;

  // Vector of composite cross sections for each material in the MaterialTable 
  // The composite cross section is composed of cross sections for each element in material
  std::vector<G4IDataSet*>* crossSections;

  std::vector<G4String> crossModel;

};

#endif











