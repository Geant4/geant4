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
// $Id: G4LEProtonInelastic.hh,v 1.12 2007-02-26 18:31:19 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy Proton Inelastic Process
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 19-Nov-1996
 // Last modified: 27-Mar-1997
 // Modified by J.L.Chuma 8-Jul-97 to implement Nucleus changes
 
#ifndef G4LEProtonInelastic_h
#define G4LEProtonInelastic_h 1
 
// Class Description
// Final state production model for Proton inelastic scattering below 20 GeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "G4InelasticInteraction.hh"
 
 class G4LEProtonInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEProtonInelastic() : G4InelasticInteraction("G4LEProtonInelastic")
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 55.*GeV );
    }
    
    ~G4LEProtonInelastic()
    { }
    
    G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack,
                                      G4Nucleus &targetNucleus );
 private:
    
    void Cascade(                               // derived from CASP
      G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
      G4int &vecLen,
      const G4HadProjectile *originalIncident,
      G4ReactionProduct &currentParticle,
      G4ReactionProduct &targetParticle,
      G4bool &incidentHasChanged, 
      G4bool &targetHasChanged,
      G4bool &quasiElastic );

    void SlowProton(
      const G4HadProjectile *originalIncident,
      G4Nucleus &targetNucleus );
 };
 
#endif
 
