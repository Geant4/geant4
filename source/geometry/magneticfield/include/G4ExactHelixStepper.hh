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
// $Id: G4ExactHelixStepper.hh 66356 2012-12-18 09:02:32Z gcosmo $ 
//
//
// class G4ExactHelixStepper
//       -------------------
// Class description:
//
// Concrete class for particle motion in constant magnetic field.

// History:
// - 28.Jan.05  J.Apostolakis   Creation of new concrete class
// --------------------------------------------------------------------

#ifndef G4ExactHelixStepper_hh
#define G4ExactHelixStepper_hh

#include "G4Types.hh"
#include "G4ThreeVector.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4MagHelicalStepper.hh"
#include "G4Mag_EqRhs.hh"

class G4ExactHelixStepper : public G4MagHelicalStepper
{
  public:  // with description

    G4ExactHelixStepper(G4Mag_EqRhs *EqRhs);
    ~G4ExactHelixStepper();
  
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[]  );
      // Step 'integration' for step size 'h'
      // Provides helix starting at y[0 to 6]
      // Outputs yout[] and ZERO estimated error yerr[]=0.
  
    void DumbStepper( const G4double y[],
                            G4ThreeVector   Bfld,
                            G4double  h,
                            G4double yout[] );
      // Performs a 'dump' Step without error calculation.
  
    G4double DistChord() const;
      // Estimate maximum distance of curved solution and chord ... 

    virtual G4int IntegratorOrder() const;

  private:

    G4ExactHelixStepper(const G4ExactHelixStepper&);
    G4ExactHelixStepper& operator=(const G4ExactHelixStepper&);
      // Private copy constructor and assignment operator.
   
  private:

    G4ThreeVector    fBfieldValue;
      //  Initial value of field at last step
    G4Mag_EqRhs*  fPtrMagEqOfMot;
};



#endif  /* G4ExactHelixStepper_hh */
