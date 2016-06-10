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
// $Id: G4ExplicitEuler.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4ExplicitEuler
//
// Class description:
//
// Explicit Euler stepper for magnetic field: x_1 = x_0 + h * dx_0.
// The most simple approach for solving linear differential equations.
// Take the current derivative and add it to the current position.

// History:
// - Created. W.Wander <wwc@mit.edu>, 12/09/97
// - Adjusted for new Equation classes, J.Apostolakis 5/11/03
// -------------------------------------------------------------------

#ifndef G4EXPLICITEULER_HH
#define G4EXPLICITEULER_HH

#include "G4MagErrorStepper.hh"

class G4ExplicitEuler : public G4MagErrorStepper
{

  public:  // with description

    G4ExplicitEuler(G4EquationOfMotion* EqRhs, G4int numberOfVariables = 6) ;
   ~G4ExplicitEuler();

    void  DumbStepper(  const G4double y[],
                        const G4double dydx[],
                              G4double h,
                              G4double yout[]);

  public:  // without description

    G4int IntegratorOrder() const { return 1; }

};

#endif /* G4EXPLICITEULER_HH */
