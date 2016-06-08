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
  ~G4NeutronIsotopeProduction();

  G4IsoResult * GetIsotope(const G4Track & aTrack, const G4Nucleus & aNucleus);

  private:
    
  G4ElementIsoCrossSections<G4NeutronIsoIsoCrossSections> ** theData;
  G4int numberOfElements;
};

#endif
