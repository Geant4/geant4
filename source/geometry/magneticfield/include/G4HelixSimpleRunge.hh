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
// $Id: G4HelixSimpleRunge.hh,v 1.6 2003/10/31 14:35:51 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
// class G4HelixSimpleRunge
//
// Class description:
//
//  Helix Simple Runge-Kutta stepper for magnetic field:
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
//  Second order solver.
//  Take the derivative at a position to be assumed at the middle of the
//  Step and add it to the current position.

// W. Wander <wwc@mit.edu> 03/12/98
// -------------------------------------------------------------------

#ifndef G4HELIXSIMPLERUNGE_HH
#define G4HELIXSIMPLERUNGE_HH
#include "G4MagHelicalStepper.hh"

class G4HelixSimpleRunge : public G4MagHelicalStepper
{
  public:  // with description

    G4HelixSimpleRunge(G4Mag_EqRhs *EqRhs)
      : G4MagHelicalStepper(EqRhs){}

    ~G4HelixSimpleRunge(){}
  
    void  DumbStepper(  const G4double y[],
                              G4ThreeVector   Bfld,
                              G4double        h,
                              G4double        yout[]);

  public:  // without description
  
    G4int IntegratorOrder() const { return 2; }
};

#endif /* G4HELIXSIMPLERUNGE_HH */
