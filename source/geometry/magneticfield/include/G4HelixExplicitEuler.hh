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
//
// $Id: G4HelixExplicitEuler.hh,v 1.6 2001-07-11 09:59:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4HelixExplicitEuler
//
// Class description:
//
// Helix Explicit Euler: x_1 = x_0 + helix(h)
// with helix(h) being a helix piece of length h.
// A simple approach for solving linear differential equations.
// Take the current derivative and add it to the current position.

// History:
// - Created. W.Wander <wwc@mit.edu>, 12/09/97

#ifndef G4HELIXEXPLICITEULER_HH
#define G4HELIXEXPLICITEULER_HH

#include "G4MagHelicalStepper.hh"

class G4HelixExplicitEuler : public G4MagHelicalStepper
{

  public:

    G4HelixExplicitEuler(G4Mag_EqRhs *EqRhs)
      : G4MagHelicalStepper(EqRhs) {;}
 
    ~G4HelixExplicitEuler() {;}
  
    void DumbStepper( const G4double y[],
		      G4ThreeVector  Bfld,
		      G4double       h,
		      G4double       yout[]);

  public:  // without description

    // DELETED  RightHandSide( ) !!!!  
    // Replaced by MagFieldEvaluate( const G4double y[], G4double B[] )   
    // in G4HelicalStepper
  
    G4int IntegratorOrder() const { return 1; }
};

#endif /* G4EXPLICITEULER_HH */
