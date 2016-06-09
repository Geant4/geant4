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
#ifndef G4ProtonIsotopeProduction_h
#define G4ProtonIsotopeProduction_h

#include "globals.hh"
#include "G4VIsotopeProduction.hh"
#include "G4ElementIsoCrossSections.hh"
#include "G4ProtonIsoIsoCrossSections.hh"
#include "Randomize.hh"

// Class Description
// Proton-induced isotope production model for energies below 100 MeV. 
// It runs in parasitic mode to the LEP inelastic models.  In your
// physics list, this class must be registered with the
// G4LEProtonInelastic model.  This class is a prototype only as the
// corresponding isotope production data does not exist in the G4NDL
// data library.
// Class Description - End

class G4ProtonIsotopeProduction : public G4VIsotopeProduction
{
  public: 
    G4ProtonIsotopeProduction() {
      G4cout << "WARNING: G4ProtonIsotopeProduction is deprecated and will be removed with Geant4 version 10" 
             << G4endl;
    }
    ~G4ProtonIsotopeProduction();

    G4IsoResult* GetIsotope(const G4HadProjectile* aTrack,
                            const G4Nucleus& aNucleus);

  private:
    G4ElementIsoCrossSections<G4ProtonIsoIsoCrossSections>** theData;
    G4int numberOfElements;
};

#endif
