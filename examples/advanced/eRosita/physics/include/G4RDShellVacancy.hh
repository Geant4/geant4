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
//  21 Sept 2001  Elena Guardincerri   Created
//  25 Mar  2002  V.Ivanchenko         Change AverageNOfIonisations int->double
//  12 Apr  2003  V.Ivanchenko         Migrade to cut per region
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Fluorescence in along step do it: determinates the number of ionizations

// -------------------------------------------------------------------

#ifndef G4RDSHELLVACANCY_HH
#define G4RDSHELLVACANCY_HH 1

#include "globals.hh"
#include <vector>

class G4RDVEMDataSet;
class G4MaterialCutsCouple;
class G4Element;
class G4RDShellVacancy
{
public:

  G4RDShellVacancy();

  ~G4RDShellVacancy();

  std::vector<G4int> GenerateNumberOfIonisations(const G4MaterialCutsCouple* couple,
						   G4double incidentEnergy,
						   G4double eLoss) const;

  void AddXsiTable(G4RDVEMDataSet* set);

private:

  G4double AverageNOfIonisations(const G4MaterialCutsCouple* couple,
			               G4int index, 
			               G4double energy, 
			               G4double eLoss) const;
  
  std::vector<G4RDVEMDataSet*> xsis;
};

#endif




