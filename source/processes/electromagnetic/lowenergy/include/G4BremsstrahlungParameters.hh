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
// $Id: G4BremsstrahlungParameters.hh,v 1.3 2001-10-09 11:23:24 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//                          Values of the parameters from A. Forti's fit
// 12.09.01 V.Ivanchenko    Add activeZ and paramA
// 25.09.01 V.Ivanchenko    Add parameter C and change interface to B
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Set of parameters for LowEnergyBremsstrahlung
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4BREMSSTRAHLUNGPARAMETERS_HH
#define G4BREMSSTRAHLUNGPARAMETERS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "g4std/map"
#include "g4std/vector"

class G4VEMDataSet;
class G4VDataSetAlgorithm;

class G4BremsstrahlungParameters {
 
public:

  G4BremsstrahlungParameters(G4int minZ = 1, G4int maxZ = 99);

  ~G4BremsstrahlungParameters();
 
  G4double ParameterA(G4int Z, G4double energy) const;

  G4double ParameterB(G4int Z, G4double energy) const;

  G4double ParameterC(G4int index) const;
  
  void PrintData() const;

private:

  // hide assignment operator 
  G4BremsstrahlungParameters(const G4BremsstrahlungParameters&);
  G4BremsstrahlungParameters & operator=(
                             const G4BremsstrahlungParameters &right);

  void LoadData();

  // The interpolation algorithm
  G4VDataSetAlgorithm* interpolation;

  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> > paramA;

  G4DataVector paramB0;   
  G4DataVector paramB1;
  G4DataVector paramC;
  G4DataVector activeZ;
  
  G4int zMin;
  G4int zMax;

};
 
#endif
 










