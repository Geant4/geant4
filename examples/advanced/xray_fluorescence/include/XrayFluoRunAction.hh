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
// $Id: XrayFluoRunAction.hh
// GEANT4 tag $Name:  xray_fluo-V04-01-03
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
#include "g4std/map"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class XrayFluoDataSet;
class G4DataVector;

class XrayFluoRunAction : public G4UserRunAction
{
  public:
   

 XrayFluoRunAction();

 ~XrayFluoRunAction();
  const XrayFluoDataSet* GetSet();
  const XrayFluoDataSet* GetGammaSet();
  const XrayFluoDataSet* GetAlphaSet();
  const XrayFluoDataSet* GetEfficiencySet();
  G4DataVector* GetEnergies();
  G4DataVector* GetData();
  

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  
 //returns the sum of the element of data
  G4double GetDataSum();

  //calulates the maximum integer contained in energy
  //choose from dataMap and energyMap the set of values
  //corresponding to the calculated integer
  //returns a value of energy chosen randomly from this set 
  G4double GetInfData(G4double energy, G4double random);

 //calulates the minimum integer contained in energy
  //choose from dataMap and energyMap the set of values
  //corresponding to the calculated integer
  //returns a value of energy chosen randomly from this set 
  G4double GetSupData(G4double energy, G4double random);

private:

  const XrayFluoDataSet* dataSet;

  //stores the data of the incident gamma spectrum
  const XrayFluoDataSet* dataGammaSet;

  //stores the data of the incident alpha spectrum
  const XrayFluoDataSet* dataAlphaSet;

  //stores the data of the efficience of the detector
  const XrayFluoDataSet* efficiencySet;

  //stores the energy data of the proton and alpha spectra
  G4DataVector* energies;

//stores the data of the proton and alpha spectra
   G4DataVector* data;

  //stores the energy data (first column of the file) of the 
  //response function 
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyMap;
  
  //stores the values (second column of the file) of the 
  //response function 
G4std::map<G4int,G4DataVector*,G4std::less<G4int> > dataMap;
  
  //read the data for protons and alpha spectra
 void ReadData(G4double,G4String);

  //read the data for the response function if the detector
 void ReadResponse(const G4String& fileName);
  
 

};

#endif

