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
// $Id: G4ExplicitEuler.hh,v 1.6 2002-11-29 13:47:49 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4ExplicitEuler
//
// Class description:
//
// Explicit Euler: x_1 = x_0 + h * dx_0.
// The most simple approach for solving linear differential equations.
// Take the current derivative and add it to the current position.

// History:
// - Created. W.Wander <wwc@mit.edu>, 12/09/97

#ifndef G4EXPLICITEULER_HH
#define G4EXPLICITEULER_HH

#include "G4MagErrorStepper.hh"

class G4ExplicitEuler : public G4MagErrorStepper
{

  public:  // with description

    G4ExplicitEuler(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;
   ~G4ExplicitEuler();

    void  DumbStepper(  const G4double y[],
		        const G4double dydx[],
		              G4double h,
			      G4double yout[]);

  public:  // without description

    G4int IntegratorOrder() const { return 1; }

};

#endif /* G4EXPLICITEULER_HH */
