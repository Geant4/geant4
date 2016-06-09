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
// $Id: G4BremsstrahlungParameters.hh,v 1.9 2003/06/16 16:59:39 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

#ifndef G4BREMSSTRAHLUNGPARAMETERS_HH
#define G4BREMSSTRAHLUNGPARAMETERS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include <map>

class G4VEMDataSet;
class G4VDataSetAlgorithm;

class G4BremsstrahlungParameters {

public:

  G4BremsstrahlungParameters(const G4String& name, size_t num, G4int minZ = 1, G4int maxZ = 99);

  ~G4BremsstrahlungParameters();

  G4double Parameter(G4int parameterIndex, G4int Z, G4double energy) const;

  G4double ParameterC(G4int index) const;

  void PrintData() const;

private:

  // hide assignment operator
  G4BremsstrahlungParameters(const G4BremsstrahlungParameters&);
  G4BremsstrahlungParameters & operator=(const G4BremsstrahlungParameters &right);

  void LoadData(const G4String& name);

  std::map<G4int,G4VEMDataSet*,std::less<G4int> > param;

  G4DataVector paramC;
  G4DataVector activeZ;

  G4int zMin;
  G4int zMax;

  size_t length;

};
 
#endif











