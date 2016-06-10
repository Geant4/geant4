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
// $Id: G4RKG3_Stepper.hh 68055 2013-03-13 14:43:28Z gcosmo $
//
//
//
// class G4RKG3_Stepper
//
// Class description:
//
// Integrator Runga-Kutta Stepper from Geant3.
//
// History:
// - Created. J.Apostolakis, V.Grichine - 30.01.97
// -------------------------------------------------------------------

#ifndef G4RKG3_Stepper_hh
#define G4RKG3_Stepper_hh

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"

class G4Mag_EqRhs;

class G4RKG3_Stepper : public G4MagIntegratorStepper
{
  public:  // with description

    G4RKG3_Stepper(G4Mag_EqRhs *EqRhs);
      // Integrate over 6 variables only:  position & velocity.
      // Not implemented yet !

    ~G4RKG3_Stepper();

    void Stepper( const G4double yIn[],
                  const G4double dydx[],
                        G4double h,
                        G4double yOut[],
                        G4double yErr[]  );
      // The method which must be provided, even if less efficient.

    G4double  DistChord() const ;
 
    void StepNoErr( const G4double tIn[8],
                    const G4double dydx[6],
                          G4double Step,
                          G4double tOut[8],
                          G4double B[3] );
      // Integrator RK Stepper from G3 with only two field evaluation per 
      // Step. It is used in propagation initial Step by small substeps
      // after solution error and delta geometry considerations. 
      // B[3] is magnetic field which is passed from substep to substep.

    void StepWithEst( const G4double  tIn[8],
                      const G4double dydx[6],
                            G4double Step,
                            G4double tOut[8],
                            G4double& alpha2,
                            G4double& beta2,
                      const G4double B1[3],
                            G4double B2[3] );
      // Integrator for RK from G3 with evaluation of error in solution and delta
      // geometry based on naive similarity with the case of uniform magnetic field.
      // B1[3] is input  and is the first magnetic field values
      // B2[3] is output and is the final magnetic field values.

  public:  // without description

    G4int IntegratorOrder() const { return 4; }

  private:

    G4ThreeVector fyInitial,
                  fyMidPoint,
                  fyFinal;
   G4ThreeVector  fpInitial;
   G4ThreeVector  BfldIn;
   G4double       hStep;
};

#endif
