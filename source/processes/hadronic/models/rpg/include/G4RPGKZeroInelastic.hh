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
// $Id: G4RPGKZeroInelastic.hh 79697 2014-03-12 13:10:09Z gcosmo $
//
// Author: D. H. Wright
// Date:   18 June 2007
//

// Class Description:
// Final state production model for K0 inelastic scattering using the 
// re-parameterized Gheisha model.  This is a utility class only and
// should be used in a physics list.  

#ifndef G4RPGKZeroInelastic_h
#define G4RPGKZeroInelastic_h 1
 
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4RPGInelastic.hh"
 
 class G4RPGKZeroInelastic : public G4RPGInelastic
 {
 public:
    
    G4RPGKZeroInelastic() : G4RPGInelastic("G4RPGKZeroInelastic")
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*CLHEP::GeV );
    }
    
    ~G4RPGKZeroInelastic()
    { }
    
    G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void Cascade(                               // derived from CASK0
      G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
      G4int &vecLen,
      const G4HadProjectile *originalIncident,
      G4ReactionProduct &currentParticle,
      G4ReactionProduct &targetParticle,
      G4bool &incidentHasChanged, 
      G4bool &targetHasChanged,
      G4bool &quasiElastic );
    
 };
 
#endif
 
