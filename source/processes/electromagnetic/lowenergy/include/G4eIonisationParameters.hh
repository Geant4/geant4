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
// $Id: G4eIonisationParameters.hh,v 1.2 2001-09-16 14:13:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//                          Values of the parameters from A. Forti's fit
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 12.09.01 V.Ivanchenko    Add param and interpolation of parameters  
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Set of parameters for LowEnergyIonisation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4IONISATIONPARAMETERS_HH
#define G4IONISATIONPARAMETERS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "g4std/map"

class G4VDataSetAlgorithm;
class G4VEMDataSet;

class G4eIonisationParameters {
 
public:
  // To be completed

  G4eIonisationParameters(G4int minZ = 1, G4int maxZ = 99);

  ~G4eIonisationParameters();
 
  const G4DataVector& Parameters(G4int Z, G4int shellIndex, G4double e) const;

  G4double Parameter(G4int Z, G4int shellIndex, G4int parameterIndex,
                     G4double e) const;
  
  void PrintData() const;

private:

  void LoadData();

  // To be completed: data member(s) for parameters 

  G4int zMin;
  G4int zMax;

  // Parameters of the energy spectra
  G4DataVector activeZ;

  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> > param;

  size_t length;

  // The interpolation algorithm
  const G4VDataSetAlgorithm* interpolation;

};
 
#endif
 










