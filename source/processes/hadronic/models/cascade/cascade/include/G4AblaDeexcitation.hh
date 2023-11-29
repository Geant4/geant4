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
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   March 2023
//
// Class that does the nuclear de-excitation, after the Bertini cascade,
// by using the Abla model, and then transform the secondaries in Bertini
// objects.
//----------------------------------------------------------------------------
//
#ifndef G4ABLA_DEEXCITATION_HH
#define G4ABLA_DEEXCITATION_HH


#include "G4CascadeDeexciteBase.hh"
#include "globals.hh"

class G4VPreCompoundModel;


class G4AblaDeexcitation : public G4CascadeDeexciteBase {
  public:
    G4AblaDeexcitation();
    ~G4AblaDeexcitation() override;
    void setVerboseLevel( G4int verbose ) override;
    void deExcite( const G4Fragment& fragment, G4CollisionOutput& globalOutput ) override;
    G4AblaDeexcitation( const G4AblaDeexcitation& ) = delete;
    G4AblaDeexcitation& operator=( const G4AblaDeexcitation& ) = delete;
  private:
    G4VPreCompoundModel* theDeExcitation;
};

#endif	/* G4ABLA_DEEXCITATION_HH */
