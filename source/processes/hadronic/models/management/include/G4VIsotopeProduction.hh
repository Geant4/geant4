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
#ifndef G4VIsotopeProduction_h
#define G4VIsotopeProduction_h 1

#include "G4IsoResult.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"

// Class Description
// This is the class you inherit from, if you want to implement your special
// isotope production model, based on your production cross-sections.
// Registering it with the corresponding process, the particle flux from the
// transport models will be fed into the production model, to retrieve improved
// isotope production information.
// Class Description - End
class G4VIsotopeProduction
{
public:// With Description

  // This is the interface to implement for isotope production models.
  virtual G4IsoResult * GetIsotope(const G4Track & aTrack, const G4Nucleus & aNucleus) = 0;
  
public:// Without Description
  G4bool operator == (const G4VIsotopeProduction & aProd)
  {
    G4bool result = false;
    if(&aProd==this) result = true;
    return result;
  }
};

#endif
