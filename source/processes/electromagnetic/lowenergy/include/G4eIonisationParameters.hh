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
// $Id: G4eIonisationParameters.hh,v 1.1 2001-08-20 16:36:01 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//                          Values of the parameters from A. Forti's fit
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

class G4DataVector;

class G4IonisationParameters {
 
public:
  // To be completed

  G4IonisationParameters(G4int minZ = 1, G4int maxZ = 99);

  ~G4IonisationParameters();
 
  const G4DataVector& Parameters(G4int Z, G4int shellIndex) const;

  G4double Parameter(G4int Z, G4int shellIndex, G4int parameterIndex) const;
  
  void PrintData() const;

private:

  void LoadData();

  // To be completed: data member(s) for parameters 

  G4int zMin;
  G4int zMax;
};
 
#endif
 










