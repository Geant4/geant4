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
//
// $Id: G4MagHelicalStepper.hh 93806 2015-11-02 11:21:01Z gcosmo $
//
//
//
// class G4MagHelicalStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field

// It is used for a set of steppers which use the helix as a sort of
// 'first order' solution.
//   - Most obtain an error by breaking up the step in two
//   - G4ExactHelicalStepper does not provide an error estimate

// History:
// - 05.11.98  J.Apostolakis   Creation of new ABC 
// --------------------------------------------------------------------

#ifndef G4MagHelicalStepper_hh
#define G4MagHelicalStepper_hh

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

class G4MagHelicalStepper : public G4MagIntegratorStepper
{
  public:  // with description

    G4MagHelicalStepper(G4Mag_EqRhs *EqRhs);
    virtual ~G4MagHelicalStepper();
  
    virtual void Stepper( const G4double y[],    // VIRTUAL for ExactHelix - temporary
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[]  );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].
  
    virtual  void DumbStepper( const G4double y[],
                               G4ThreeVector   Bfld,
                               G4double  h,
                               G4double yout[] ) = 0;
      // Performs a 'dump' Step without error calculation.
  
    G4double DistChord()const ;
      // Estimate maximum distance of curved solution and chord ... 

  protected:  // with description

    inline void LinearStep( const G4double  yIn[],
                                  G4double  h,
                                  G4double  yHelix[]) const;
      // A linear Step in regions without magnetic field.

     void AdvanceHelix( const G4double  yIn[],
                             G4ThreeVector   Bfld,
                             G4double  h,
			G4double  yHelix[],G4double yHelix2[]=0);    // output 
      // A first order Step along a helix inside the field.

    inline void MagFieldEvaluate( const G4double y[], G4ThreeVector& Bfield );
      // Evaluate the field at a certain point.

  
   inline G4double GetInverseCurve( const G4double Momentum, const G4double Bmag );
      // Evaluate Inverse of Curvature of Track

      // Store and use the parameters of track : 
      // Radius of curve, Stepping angle, Radius of projected helix
   inline void SetAngCurve(const G4double Ang);
   inline G4double GetAngCurve()const;

   inline void SetCurve(const G4double Curve);
   inline G4double GetCurve()const;

   inline void SetRadHelix(const G4double Rad);
   inline G4double GetRadHelix()const;


  protected:  // without description

    // void MagFieldEvaluate( const G4double y[], G4double B[] )   
    //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

  private:

    G4MagHelicalStepper(const G4MagHelicalStepper&);
    G4MagHelicalStepper& operator=(const G4MagHelicalStepper&);
      // Private copy constructor and assignment operator.
 
    static const G4double fUnitConstant;   //  As in G4Mag_EqRhs.hh/cc where it is not used.
  private:
   
    G4Mag_EqRhs*  fPtrMagEqOfMot;
 
    // Data stored in order to find the chord.
      G4double fAngCurve;
      G4double frCurve;
      G4double frHelix;
    // Data stored in order to find the chord.
      G4ThreeVector yInitial, yMidPoint, yFinal;
       
    
};

#include  "G4MagHelicalStepper.icc"

#endif  /* G4MagHelicalStepper_hh */
