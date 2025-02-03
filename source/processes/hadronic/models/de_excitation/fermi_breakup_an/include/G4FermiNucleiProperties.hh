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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMINUCLEIPROPERTIES_HH
#define G4FERMINUCLEIPROPERTIES_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiFastNucleiProperties.hh"
#include "G4FermiSingleton.hh"
#include "G4FermiVNucleiProperties.hh"

// it is possible to use polymorphism here
// but it is a bottleneck and no virtual call is made
using G4FermiNucleiProperties = G4FermiSingleton<G4FermiFastNucleiProperties>;

static_assert(
  std::is_base_of_v<G4FermiVNucleiProperties,
                    std::remove_reference_t<decltype(G4FermiNucleiProperties::Instance())>>,
  "Incorrect Nuclei Properties class");

#endif  // G4FERMINUCLEIPROPERTIES_HH
