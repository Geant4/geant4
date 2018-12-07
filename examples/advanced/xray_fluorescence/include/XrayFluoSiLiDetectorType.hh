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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  19 Jun 2003  Alfonso Mantero Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoSiLidetectorType_hh
#define XrayFluoSiLidetectorType_hh 1

#include "globals.hh"
#include <map>
#include "G4DataVector.hh"
#include "XrayFluoDataSet.hh"
#include "XrayFluoVDetectorType.hh"

//class XrayFluoDataSet;
//class G4DataVector;

class XrayFluoSiLiDetectorType : public XrayFluoVDetectorType
{
public:
  
  ~XrayFluoSiLiDetectorType();
  
  //returns an unique instance of this class (singleton)
  static XrayFluoSiLiDetectorType* GetInstance();

  //gives the material (a string) of wich is made the detector
  G4String GetDetectorMaterial();
  
  // Given the energy depositon in the detector, returns the measured value,
  // according the energy definition power of the detector
  G4double ResponseFunction(G4double);
  
  //returns a random value of the energy measured according to the tabulated
  //peack just lower of the energy deposited
  G4double GetInfData(G4double,G4double,G4int);
  
  //returns a random value of the energy measured according to the tabulated
  //peack just higher of the energy deposited
  G4double GetSupData(G4double,G4double,G4int);

  //loads the data tabulated for the response and the efficiency data from file
  void LoadResponseData(G4String);
  void LoadEfficiencyData(G4String);
  
private:

  
  XrayFluoSiLiDetectorType();

  static XrayFluoSiLiDetectorType* instance;

  G4String detectorMaterial;

  //stores the data of the efficience of the detector
  const XrayFluoDataSet* efficiencySet;

  G4VDataSetAlgorithm* interpolation4;

  //stores the energy data (first column of the file) of the 
  //response function 
  std::map<G4int,G4DataVector*,std::less<G4int> > energyMap;
  
  //stores the values (second column of the file) of the 
  //response function 
  std::map<G4int,G4DataVector*,std::less<G4int> > dataMap;

};
#endif
