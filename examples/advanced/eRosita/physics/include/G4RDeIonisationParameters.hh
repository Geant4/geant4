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
//         V. Ivanchenko 
//         Values of the parameters from A. Forti's fit
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 12.09.01 V.Ivanchenko    Add param and interpolation of parametersVI  
// 10.10.2001 MGP           Revision to improve code quality and 
//                          consistency with design
// 29.11.01  V.Ivanchenko    Parametrisation is updated
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Set of parameters for LowEnergyIonisation described spectrum 
// of delta-electrons retrieved from EEDL database.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDIONISATIONPARAMETERS_HH
#define G4RDIONISATIONPARAMETERS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include <map>

class G4RDVDataSetAlgorithm;
class G4RDVEMDataSet;

class G4RDeIonisationParameters {
 
public:

  G4RDeIonisationParameters(G4int minZ = 1, G4int maxZ = 99);

  ~G4RDeIonisationParameters();
 
  G4double Parameter(G4int Z, G4int shellIndex, 
		     G4int parameterIndex, G4double e) const;

  G4double Excitation(G4int Z, G4double e) const;
  
  void PrintData() const;

private:

  // Hide copy constructor and assignment operator 
  G4RDeIonisationParameters(const G4RDeIonisationParameters&);
  G4RDeIonisationParameters & operator=(const G4RDeIonisationParameters &right);

  void LoadData();

  G4int zMin;
  G4int zMax;

  G4DataVector activeZ;

  // Parameters of the energy spectra
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> > param;
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> > excit;

  size_t length;
};
 
#endif
