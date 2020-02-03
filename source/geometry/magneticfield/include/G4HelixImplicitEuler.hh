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
// G4HelixImplicitEuler
//
// Class description:
//
//  Helix Implicit Euler stepper for magnetic field:
//        x_1 = x_0 + 1/2 * ( helix(h,t_0,x_0)
//                          + helix(h,t_0+h,x_0+helix(h,t0,x0) ) )
//  Second order solver.
//  Take the current derivative and add it to the current position.
//  Take the output and its derivative. Add the mean of both derivatives
//  to form the final output.

// Author: W.Wander <wwc@mit.edu>, 03/11/1998
// -------------------------------------------------------------------
#ifndef G4HELIXIMPLICITEULER_HH
#define G4HELIXIMPLICITEULER_HH

#include "G4MagHelicalStepper.hh"

class G4HelixImplicitEuler : public G4MagHelicalStepper
{

  public:  // with description

    G4HelixImplicitEuler(G4Mag_EqRhs *EqRhs);
    ~G4HelixImplicitEuler();
  
    void DumbStepper( const G4double y[],
                            G4ThreeVector Bfld,
                            G4double      h,
                            G4double      yout[] );
  
  public:  // without description

    inline G4int IntegratorOrder() const { return 2; }
};

#endif
