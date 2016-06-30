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
// $Id: G4VFSALIntegrationStepper.hh
//
// class G4VFSALIntegrationStepper
//
// Class description:
//     Class similar to G4VMagIntegratorStepper, for steppers which
//   estimate the value of the derivative at the projected endpoint
//   of integration - at each successful step.
//     This ability is known as 'First Same As Last' (FSAL). It
//   reduces the number of required calls to the equation's 
//   RightHandSide method, and, as such the number of calls to the 
//   (potentially expensive) field evaluation methods.
//
// Based on G4VMagIntegratorStepper
//
//  Design/first implementation: Somnath Banerjee, May-Aug 2015,
//     Work supported by the Google Summer of Code 2015.
// Supervision/improvement: John Apostolakis 2015-2016
// --------------------------------------------------------------------

#ifndef FSAL_MAGIntegrator_STEPPER
#define FSAL_MAGIntegrator_STEPPER

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"

class G4VFSALIntegrationStepper
{
  public:  // with description

     G4VFSALIntegrationStepper (G4EquationOfMotion* Equation,
                              G4int              numIntegrationVariables,
                              G4int              numStateVariables=12);
    
     virtual ~G4VFSALIntegrationStepper();
       // Constructor and destructor. No actions.

     virtual  void  Stepper(  const G4double y[],
                              const G4double dydx[],
                                    G4double h,
                                    G4double yout[],
                                    G4double yerr[],
		                            G4double lastDydx[]) = 0 ;
       // The stepper for the Runge Kutta integration.
       // The stepsize is fixed, with the Step size given by h.
       // Integrates ODE starting values y[0 to 6].
       // Outputs yout[] and its estimated error yerr[].

     virtual  G4double  DistChord() const = 0; 
       // Estimate the maximum distance of a chord from the true path
       // over the segment last integrated.

     virtual void ComputeRightHandSide( const G4double y[], G4double dydx[] ); 
       // Must compute the RightHandSide as in the method below
       // Optionally can cache the input y[] and the dydx[] values computed.
    
//     virtual G4bool isFSAL() const = 0;
//    	//Return true if the stepper uses FSAL (First Same As Last)
    
//    G4double *getLastDydx() {return 0;}

     inline void NormaliseTangentVector( G4double vec[6] );
       // Simple utility function to (re)normalise 'unit velocity' vector.

     inline void NormalisePolarizationVector( G4double vec[12] );
       // Simple utility function to (re)normalise 'unit spin' vector.

      void RightHandSide( const double y[], double dydx[] );
       // Utility method to supply the standard Evaluation of the
       // Right Hand side of the associated equation.


     inline G4int  GetNumberOfVariables() const;
       // Get the number of variables that the stepper will integrate over.

     // void   SetNumberOfVariables(G4int newNo);  // Dangerous & obsolete ...

     inline G4int  GetNumberOfStateVariables() const;
       // Get the number of variables of state variables (>= above, integration)

     virtual G4int IntegratorOrder() const = 0;
       // Returns the order of the integrator
       // i.e. its error behaviour is of the order O(h^order).

     inline G4EquationOfMotion *GetEquationOfMotion(); 
       // As some steppers (eg RKG3) require other methods of Eq_Rhs
       // this function allows for access to them.
     inline void SetEquationOfMotion(G4EquationOfMotion* newEquation); 

    //--- --- For DEBUG --- ---
    inline G4int GetfNoRHSCalls(){
        return fNoRHSCalls;
    }
    void increasefNORHSCalls();
    
    inline void ResetfNORHSCalls(){
        fNoRHSCalls = 0;
    }
    //--- --- ///////// --- ---
  private:
  
     G4VFSALIntegrationStepper(const G4VFSALIntegrationStepper&);
     G4VFSALIntegrationStepper& operator=(const G4VFSALIntegrationStepper&);
       // Private copy constructor and assignment operator.

  private:

     G4EquationOfMotion *fEquation_Rhs;
     const G4int  fNoIntegrationVariables;  // Number of Variables in integration
     const G4int  fNoStateVariables;        // Number required for FieldTrack
    
    
    //--- --- For DEBUG --- ---
    G4int fNoRHSCalls;
    //--- --- ///////// --- ---
    
     // const G4int  fNumberOfVariables;

};

#include  "G4VFSALIntegrationStepper.icc"

#endif  /* G4VFSALIntegrationStepper */
