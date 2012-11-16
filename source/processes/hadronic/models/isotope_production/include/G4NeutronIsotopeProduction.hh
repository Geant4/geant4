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
// Neutron-induced isotope production model for E_n < 100 MeV. 
// It runs in parasitic mode to the LEP inelastic models.
// In your physics list, this class must be registered with 
// the G4LENeutronInelastic model.  The corresponding isotope
// production data from the G4NDL neutron data library is also
// required.
// Class Description - End

class G4NeutronIsotopeProduction : public G4VIsotopeProduction
{
  public:
    G4NeutronIsotopeProduction();
    G4NeutronIsotopeProduction(const G4NeutronIsotopeProduction& nip);
    virtual ~G4NeutronIsotopeProduction();

    G4IsoResult* GetIsotope(const G4HadProjectile* aTrack,
                            const G4Nucleus& aNucleus);

  private:
    G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections>** theData;
    G4int numberOfElements;
};

#endif
