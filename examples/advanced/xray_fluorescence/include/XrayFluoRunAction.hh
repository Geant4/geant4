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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------


#ifndef XrayFluoRunAction_h
#define XrayFluoRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class XrayFluoDataSet;
class G4DataVector;

class XrayFluoRunAction : public G4UserRunAction
{

public:
  XrayFluoRunAction();
  
  ~XrayFluoRunAction();
  const XrayFluoDataSet* GetSet() const;
  const XrayFluoDataSet* GetGammaSet() const;
  const XrayFluoDataSet* GetAlphaSet() const;
  
  //  const XrayFluoDataSet* GetEfficiencySet();
  G4DataVector* GetEnergies() const;
  G4DataVector* GetData() const;
  
  
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  
  //returns the sum of the element of data
  G4double GetDataSum() const;
  
  //calulates the maximum integer contained in energy
  //choose from dataMap and energyMap the set of values
  //corresponding to the calculated integer
  //returns a value of energy chosen randomly from this set 
  //  G4double GetInfData(G4double energy, G4double random, G4int posIndex);
  
  //calulates the minimum integer contained in energy
  //choose from dataMap and energyMap the set of values
  //corresponding to the calculated integer
  //returns a value of energy chosen randomly from this set 
  //  G4double GetSupData(G4double energy, G4double random, G4int posIndex);
  
private:

  void Initialise();
  G4bool isInitialized;

  const XrayFluoDataSet* dataSet;
  
  //stores the data of the incident gamma spectrum
  const XrayFluoDataSet* dataGammaSet;

  //stores the data of the incident alpha spectrum
  const XrayFluoDataSet* dataAlphaSet;  

  //stores the energy data of the proton and alpha spectra
  G4DataVector* energies;

  //stores the data of the proton and alpha spectra
   G4DataVector* data;
  
  //read the data for protons and alpha spectra
  void ReadData(G4double,G4String);

  //read the data for the response function if the detector
  //void ReadResponse(const G4String& fileName);
  
 

};

#endif

