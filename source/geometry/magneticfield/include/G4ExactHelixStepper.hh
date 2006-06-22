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
// $Id: G4ExactHelixStepper.hh,v 1.3 2006-06-22 10:20:37 japost Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
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

    // static const G4double fUnitConstant;   // ???????????
            //  As in G4Mag_EqRhs.hh/cc where it is not used.
  private:
    G4ThreeVector    fBfieldValue;   //  Initial value of field at last step
    G4ThreeVector    yInitialEHS,  yFinalEHS;  
    G4double         fLastStepSize;  // Length of last step
    G4double         fYInSav[7];     // Starting state of  x, p, ...
     // Values saved for calculating mid-point for chord

    // G4Mag_EqRhs*  fPtrMagEqOfMot;
};

// #include  "G4ExactHelixStepper.icc"

#endif  /* G4ExactHelixStepper_hh */
