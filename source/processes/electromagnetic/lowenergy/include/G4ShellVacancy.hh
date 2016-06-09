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
// $Id: G4ShellVacancy.hh
// GEANT4 tag $Name: 
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

#ifndef G4SHELLVACANCY_HH
#define G4SHELLVACANCY_HH 1

#include "globals.hh"
#include <vector>

class G4VEMDataSet;
class G4MaterialCutsCouple;
class G4Element;
class G4ShellVacancy
{
public:

  G4ShellVacancy();

  ~G4ShellVacancy();

  std::vector<G4int> GenerateNumberOfIonisations(const G4MaterialCutsCouple* couple,
						   G4double incidentEnergy,
						   G4double eLoss) const;

  void AddXsiTable(G4VEMDataSet* set);

private:

  G4double AverageNOfIonisations(const G4MaterialCutsCouple* couple,
			               G4int index, 
			               G4double energy, 
			               G4double eLoss) const;
  
  std::vector<G4VEMDataSet*> xsis;
};

#endif




