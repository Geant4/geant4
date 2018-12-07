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
// class G4NystromRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// Notes: 1) field must be time-independent.
//        2) time is not integrated
// 
// History:
// - Created: I.Gavrilenko   15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#ifndef G4NYSTROMRK4_HH
#define G4NYSTROMRK4_HH

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4CachedMagneticField.hh"
#include "G4ThreeVector.hh"

#include <memory>

class G4NystromRK4 : public G4MagIntegratorStepper
{
public: 
    // Can be used only for Magnetic Fields - and for 6 variables (x,p)
    G4NystromRK4(G4Mag_EqRhs* EquationMotion, 
                 G4double distanceConstField = 0.0); 

    // Single call for integration result and error
    // - Provides Error via analytical method
    virtual void Stepper(const G4double y[],
                         const G4double dydx[],
                         G4double hstep,
                         G4double yOut[],
                         G4double yError[]) override;

    void SetDistanceForConstantField(G4double length); 
    G4double GetDistanceForConstantField() const; 
   
    virtual G4int IntegratorOrder() const override { return 4; }
    virtual G4double DistChord() const override; 
  
private:
    inline void GetFieldValue(const G4double point[4], G4double field[3]);
    inline G4double GetFCof();

    G4CachedMagneticField* GetField();
    const G4CachedMagneticField* GetField() const;

    G4double fMomentum;
    G4double fMomentum2;
    G4double fInverseMomentum;
    G4double fCoefficient;
    G4ThreeVector fInitialPoint;
    G4ThreeVector fMidPoint;
    G4ThreeVector fEndPoint;

    std::unique_ptr<G4CachedMagneticField> fCachedField;
};

#include "G4NystromRK4.icc"

#endif
