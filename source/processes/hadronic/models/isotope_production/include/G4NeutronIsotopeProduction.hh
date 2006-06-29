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
#ifndef G4NeutronIsotopeProduction_h
#define G4NeutronIsotopeProduction_h

#include "globals.hh"
#include "G4VIsotopeProduction.hh"
#include "G4ElementIsoCrossSections.hh"
#include "G4NeutronIsoIsoCrossSections.hh"
#include "Randomize.hh"

// Class Description
// Isotope production model for neutron induced production (E_n<100MeV); 
// Runs in parasitic mode to the transport models.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process. You will also need the corresponding
// isotope production data from the neutron data library.
// Class Description - End

class G4NeutronIsotopeProduction : public G4VIsotopeProduction
{
  public:
  
  G4NeutronIsotopeProduction();
  virtual ~G4NeutronIsotopeProduction();

  G4IsoResult * GetIsotope(const G4Track & aTrack, const G4Nucleus & aNucleus);

  private:
    
  G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections> ** theData;
  G4int numberOfElements;
};

#endif
