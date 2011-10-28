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
// $Id: G4SeltzerBergerModel.cc,v 1.18 2010-11-04 17:30:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SeltzerBergerModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 12.08.2008
//
// Modifications:
//
// 13.11.08    add SetLPMflag and SetLPMconstant methods
// 13.11.08    change default LPMconstant value
// 13.10.10    add angular distributon interface (VI)
//
// Main References:
//  Y.-S.Tsai, Rev. Mod. Phys. 46 (1974) 815; Rev. Mod. Phys. 49 (1977) 421. 
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media, Wiley, 1972.
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4SeltzerBergerModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4ModifiedTsai.hh"

#include "G4Physics2DVector.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4SeltzerBergerModel::G4SeltzerBergerModel(const G4ParticleDefinition* p,
					   const G4String& name)
  : G4eBremsstrahlungRelModel(p,name)
{
  SetLowEnergyLimit(0.0);
  SetLPMFlag(false);
  dataSB.resize(101,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SeltzerBergerModel::~G4SeltzerBergerModel()
{
  for(size_t i=0; i<101; ++i) { delete dataSB[i]; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SeltzerBergerModel::Initialise(const G4ParticleDefinition* p,
				      const G4DataVector& cuts)
{
  // check environment variable
  // Build the complete string identifying the file with the data set
  char* path = getenv("G4BREM");

  // Access to elements
  const G4ElementTable* theElmTable = G4Element::GetElementTable();
  size_t numOfElm = G4Element::GetNumberOfElements();
  if(numOfElm > 0) {
    for(size_t i=0; i<numOfElm; ++i) {
      G4int Z = G4int(((*theElmTable)[i])->GetZ());
      if(Z < 1)        { Z = 1; }
      else if(Z > 100) { Z = 100; }
      G4cout << "Z= " << Z << G4endl;
      // Initialisation
      if(!dataSB[Z]) { ReadData(Z, path); }
    }
  }

  G4eBremsstrahlungRelModel::Initialise(p, cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SeltzerBergerModel::ReadData(size_t Z, const char* path)
{
  //  G4cout << "ReadData Z= " << Z << G4endl;
  // G4cout << "Status for Z= " << dataSB[Z] << G4endl;
  //if(path) { G4cout << path << G4endl; }
  if(dataSB[Z]) { return; }
  const char* datadir = path;
  if(!datadir) {
    datadir = getenv("G4BREM");
  }
  G4Physics2DVector* v = new G4Physics2DVector();
  std::ostringstream ost;
  ost << datadir << "/br" << Z;
  std::ifstream fin(ost.str().c_str());
  G4cout << "G4SeltzerBergerModel read from <" << ost.str().c_str() 
	 << ">" << G4endl;
  if(v->Retrieve(fin)) { dataSB[Z] = v; }
  // G4cout << dataSB[Z] << G4endl;
  std::ofstream* fout = new std::ofstream("test");
  v->Store(*fout);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4SeltzerBergerModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
{

  if(gammaEnergy <= 0.0 || kinEnergy <= 0.0) { return 0.0; }
  G4double x = gammaEnergy/kinEnergy;
  G4double y = log(kinEnergy/MeV);
  G4int Z = G4int(currentZ);
  //G4cout << "G4SeltzerBergerModel::ComputeDXSectionPerAtom Z= " << Z
  //	 << " x= " << x << " y= " << y << " " << dataSB[Z] << G4endl;
  if(!dataSB[Z]) { ReadData(Z); }
  G4double cross = dataSB[Z]->Value(x,y)*totalEnergy*totalEnergy
    *millibarn/(bremFactor*kinEnergy*(kinEnergy + 2*electron_mass_c2));

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


