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
// $Id: G4HelixHeum.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4HelixHeum
//
// Class description:
//
// Simple Heum stepper for magnetic field:
//     x_1 = x_0 +
//           h * 1/4 * dx(t0,x0)  +
//               3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
// Third order solver.

// History:
// - Created. W.Wander <wwc@mit.edu>, 03/11/98
// -------------------------------------------------------------------

#ifndef G4HELIXHEUM_HH
#define G4HELIXHEUM_HH

#include "G4MagHelicalStepper.hh"

class G4HelixHeum : public G4MagHelicalStepper
{

  public:  // with description

    G4HelixHeum(G4Mag_EqRhs *EqRhs)
      : G4MagHelicalStepper(EqRhs) {}

    ~G4HelixHeum() {}
  
    void DumbStepper( const G4double y[],
                            G4ThreeVector  Bfld,
                            G4double       h,
                            G4double       yout[]);

  public: // without description
  
    G4int IntegratorOrder() const { return 2; }
};

#endif /* G4HELIXHEUM_HH */
