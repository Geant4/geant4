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
// $Id: G4ConstRK4.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// Class G4ConstRK4
//
// class description:
//
// G4ConstRK4 performs the integration of one step with error calculation
// in constant magnetic field. The integration method is the same as in
// ClassicalRK4. The field value is assumed constant for the step.
// This field evaluation is called only once per step.
// G4ConstRK4 can be used only for magnetic fields.

// History:
// - 18.09.2008 - J.Apostolakis, T.Nikitina - Created
// -------------------------------------------------------------------

#ifndef G4CONSTRK4_HH
#define G4CONSTRK4_HH

#include "G4MagErrorStepper.hh"
#include "G4EquationOfMotion.hh"
#include "G4Mag_EqRhs.hh"

class G4ConstRK4 : public G4MagErrorStepper 
{
   public:  // with description

    G4ConstRK4(G4Mag_EqRhs *EquationMotion, G4int numberOfStateVariables=8);
    ~G4ConstRK4();

     void Stepper( const G4double y[],
                   const G4double dydx[],
                         G4double h,
                         G4double yout[],
                         G4double yerr[]  );
     void DumbStepper( const G4double  yIn[],
                       const G4double  dydx[],
                             G4double  h,
                             G4double  yOut[] ) ;
     G4double DistChord() const;   
 
     inline void  RightHandSideConst(const  G4double y[],
                                            G4double dydx[] ) const;

     inline void  GetConstField(const G4double y[],G4double Field[]);

   public:  // without description

     G4int IntegratorOrder() const { return 4; }

   private:

     G4ConstRK4(const G4ConstRK4&);
     G4ConstRK4& operator=(const G4ConstRK4&);
       // Private copy constructor and assignment operator.

   private:

     G4ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
     // Data stored in order to find the chord
     G4double *dydxm, *dydxt, *yt; // scratch space - not state 
     G4double *yInitial, *yMiddle, *dydxMid, *yOneStep;
     G4Mag_EqRhs *fEq;
     G4double Field[3];
};

// Inline methods

inline void G4ConstRK4:: RightHandSideConst(const G4double y[],
                                                  G4double dydx[] ) const
{
  
  G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );
    
  G4double cof =fEq->FCof()*inv_momentum_magnitude;

  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V
 
  dydx[3] = cof*(y[4]*Field[2] - y[5]*Field[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*Field[0] - y[3]*Field[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*Field[1] - y[4]*Field[0]) ;   // Az = a*(Vx*By - Vy*Bx)
}

inline void G4ConstRK4::GetConstField(const G4double y[],G4double B[])
{
  G4double  PositionAndTime[4];

  PositionAndTime[0] = y[0];
  PositionAndTime[1] = y[1];
  PositionAndTime[2] = y[2];
  // Global Time
  PositionAndTime[3] = y[7];  
  fEq -> GetFieldValue(PositionAndTime, B) ;
}

#endif  // G4CONSTRK4_HH
