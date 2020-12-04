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
//---------------------------------------------------------------------------
//
// ClassName:    G4CRCoalescence   ("CR" stands for "Cosmic Ray")
//
// Author:       2020 Alberto Ribon , based on code written by
//               Diego Mauricio Gomez Coral for the GAPS Collaboration
//
// Description:  This class can be optionally used in the method:
//
//                 G4TheoFSGenerator::ApplyYourself
//
//               to coalesce pairs of proton-neutron and antiproton-antineutron
//               into deuterons and antideuterons, respectively, from the list
//               of secondaries produced by a string model.
//               This class can be useful in particular for Cosmic Ray (CR)
//               applications.
//               By default, this class is not used.
//               However, it can be enabled via the UI command:
//
//                 /process/had/enableCRCoalescence true
//
//               It is assumed that the candidate proton-neutron and
//               antiproton-antideuteron pairs originate from the same
//               spatial position, so the condition for coalescence takes
//               into account only their closeness in momentum space.
//
//               This class is based entirely on code written by
//               Diego Mauricio Gomez Coral for the GAPS Collaboration.
//               The main application of this work is for cosmic ray physics.
//
//               Notes:
//               -  In its current version, coalescence can occur only for
//                  proton projectile (because the coalescence parameters
//                  for deuteron and antideuteron are set to non-null values
//                  only for the case of proton projectile).
//               -  This class is not meant be used for secondaries produces
//                  by intranuclear cascade models - such as BERT, BIC and
//                  INCL - which should have already a coalescence phase.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4CRCoalescence_h
#define G4CRCoalescence_h 1

#include "G4ReactionProductVector.hh"
#include "G4HadProjectile.hh"
#include "G4HadronicInteraction.hh"

class G4CRCoalescence : public G4HadronicInteraction {
  public:
  
    explicit G4CRCoalescence();
    ~G4CRCoalescence() override;
    G4CRCoalescence( const G4CRCoalescence &right ) = delete;
    const G4CRCoalescence & operator=( const G4CRCoalescence &right ) = delete;
    G4bool operator==( const G4CRCoalescence &right ) const = delete;
    G4bool operator!=( const G4CRCoalescence &right ) const = delete;

    // Set the parameter used in the coalescence condition
    void SetP0Coalescence( const G4HadProjectile &thePrimary, G4String /* model */ );

    // Main method: form deuterons and antideuterons by coalescence of, respectively,
    //              proton-neutron and antiproton-antineutron pairs with close momenta
    void GenerateDeuterons( G4ReactionProductVector* result );

  private:
 
    // Utility methods
    void PushDeuteron( const G4ThreeVector &p1, const G4ThreeVector &p2, G4int charge,
		       G4ReactionProductVector* result );
    G4int FindPartner( const G4ThreeVector &p1, G4double m1,
  		       std::vector< std::pair< G4int, G4ThreeVector > > &neutron,
		       G4double m2, G4int charge );
    G4bool Coalescence( const G4ThreeVector &p1, G4double m1,
			const G4ThreeVector &p2, G4double m2, G4int charge );
    G4bool Coalescence( G4double p1x, G4double p1y, G4double p1z, G4double m1,
                        G4double p2x, G4double p2y, G4double p2z, G4double m2, G4int charge );
    G4double GetPcm( const G4ThreeVector& p1, G4double m1,
		     const G4ThreeVector& p2, G4double m2 );
    G4double GetPcm( G4double p1x, G4double p1y, G4double p1z, G4double m1,
		     G4double p2x, G4double p2y, G4double p2z, G4double m2 );
    G4double GetS( G4double p1x, G4double p1y, G4double p1z, G4double m1,
		   G4double p2x, G4double p2y, G4double p2z, G4double m2 );
  
    G4double fP0_d;     // Coalescence parameter for deuterons
    G4double fP0_dbar;  // Coalescence parameter for antideuterons 
};

#endif
