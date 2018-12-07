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
//         V. Ivanchenko (Vladimir.Ivantchenko@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//                          Values of the parameters from A. Forti's fit
// 12.09.01 V.Ivanchenko    Add activeZ and paramA
// 25.09.01 V.Ivanchenko    Add parameter C and change interface to B
// 29.11.01 V.Ivanchenko    Parametrisation is updated
// 21.02.03 V.Ivanchenko    Number of parameters is defined in the constructor
// 28.02.03 V.Ivanchenko    Filename is defined in the constructor
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Load and access to parameters for LowEnergyBremsstrahlung from EEDL
// database. Parametrisation is described in Physics Reference Manual
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDBREMSSTRAHLUNGPARAMETERS_HH
#define G4RDBREMSSTRAHLUNGPARAMETERS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include <map>

class G4RDVEMDataSet;
class G4RDVDataSetAlgorithm;

class G4RDBremsstrahlungParameters {

public:

  G4RDBremsstrahlungParameters(const G4String& name, size_t num, G4int minZ = 1, G4int maxZ = 99);

  ~G4RDBremsstrahlungParameters();

  G4double Parameter(G4int parameterIndex, G4int Z, G4double energy) const;

  G4double ParameterC(G4int index) const;

  void PrintData() const;

private:

  // hide assignment operator
  G4RDBremsstrahlungParameters(const G4RDBremsstrahlungParameters&);
  G4RDBremsstrahlungParameters & operator=(const G4RDBremsstrahlungParameters &right);

  void LoadData(const G4String& name);

  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> > param;

  G4DataVector paramC;
  G4DataVector activeZ;

  G4int zMin;
  G4int zMax;

  size_t length;

};
 
#endif











