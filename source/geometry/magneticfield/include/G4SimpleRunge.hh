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
// $Id: G4SimpleRunge.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4SimpleRunge
//
// Class description:
//
// Simple Runge:
//
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
// Second order solver.
// Takes the derivative at a position to be assumed at the middle of the
// Step and adds it to the current position.

// History:
// - Created. W.Wander <wwc@mit.edu>, 12/09/97
// -------------------------------------------------------------------

#ifndef G4SIMPLERUNGE_HH
#define G4SIMPLERUNGE_HH

#include "G4MagErrorStepper.hh"

class G4SimpleRunge : public G4MagErrorStepper
{

  public:  // with description

    G4SimpleRunge(G4EquationOfMotion *EquationRhs, G4int numberOfVariables = 6) ;
   ~G4SimpleRunge();
      // Constructor and destructor.

    void DumbStepper( const G4double y[],
                      const G4double dydx[],
                            G4double h,
                            G4double yout[]);

  public:  // without description
  
    G4int IntegratorOrder() const { return 2; }

  private:

    G4int fNumberOfVariables ;

    G4double* dydxTemp;
    G4double* yTemp;
      // scratch space    
};

#endif /* G4SIMPLERUNGE_HH */
