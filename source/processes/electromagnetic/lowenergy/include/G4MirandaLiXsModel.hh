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
// History:
// -----------
//  30 Mar 2014   A.M., S.I. - 1st implementation
// 
// Class description
// ----------------
//  Computation of L shell experimental ionisation cross sections for protons
//  Atomic Data and Nuclear Data Tables 100 (2014) 651–780
// ---------------------------------------------------------------------------------------

#ifndef G4MirandaLiXsModel_HH
#define G4MirandaLiXsModel_HH 1

#include "globals.hh"
#include <map>

class G4VDataSetAlgorithm;
class G4VEMDataSet;

class G4MirandaLiXsModel
{
public:

  G4MirandaLiXsModel();

  virtual ~G4MirandaLiXsModel();
			     
  G4double CalculateL1CrossSection (G4int zTarget, G4double mass, G4double energyIncident);
  G4double CalculateL2CrossSection (G4int zTarget, G4double mass, G4double energyIncident);
  G4double CalculateL3CrossSection (G4int zTarget, G4double mass, G4double energyIncident);				     

private:

  G4MirandaLiXsModel(const G4MirandaLiXsModel&);
  G4MirandaLiXsModel & operator = (const G4MirandaLiXsModel &right);

  G4VDataSetAlgorithm* interpolation;

  std::map< G4int , G4VEMDataSet* > protonL1DataSetMap;
  std::map< G4int , G4VEMDataSet* > protonL2DataSetMap;
  std::map< G4int , G4VEMDataSet* > protonL3DataSetMap;

};

#endif
