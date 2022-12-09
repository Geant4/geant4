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
// G4BorisScheme
//
// Class description:
//
// Implementation of the Boris algorithm for advancing 
// charged particles in an electromagnetic field.

// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun
// --------------------------------------------------------------------
#ifndef G4BORIS_SCHEME_HH
#define G4BORIS_SCHEME_HH

class G4EquationOfMotion;

#include "G4Types.hh"


// class G4EqMagElectricField;

// #include "G4FieldTrack.hh"

#include <CLHEP/Units/PhysicalConstants.h>

class G4BorisScheme
{
  public:

    G4BorisScheme() = default;
    G4BorisScheme( // G4EqMagElectricField
       G4EquationOfMotion* equation,
                        G4int nvar = 6);
   ~G4BorisScheme() = default;

    void DoStep( G4double restMass, G4double charge, const G4double yIn[], 
                 G4double yOut[], G4double hstep) const;

  protected:
    // Used to implement the 'DoStep' method above
    void UpdatePosition(const G4double restMass, const G4double charge, const G4double yIn[],
                 G4double yOut[], G4double hstep) const;

    void UpdateVelocity(const G4double restMass, const G4double charge, const G4double yIn[],
                 G4double yOut[], G4double hstep) const;

  public: 
    // - Methods using the Boris Scheme Stepping to estimate integration error
    void StepWithErrorEstimate(const G4double yIn[], G4double restMass, G4double charge, G4double hstep,
                               G4double yOut[], G4double yErr[]) const;
      // Use two half-steps (comparing to a full step) to obtain output and error estimate

   void StepWithMidAndErrorEstimate(const G4double yIn[], G4double restMass, G4double charge, G4double hstep,
                               G4double yMid[], G4double yOut[], G4double yErr[]) const;
      // Same, and also return mid-point evaluation

    // Auxiliary method
    inline G4EquationOfMotion* GetEquationOfMotion();
    // inline void SetEquationOfMotion(G4EquationOfMotion* equation);  // Un-needed, dangerous

    inline G4int GetNumberOfVariables() const;
   
  private:

    void copy(G4double dst[], const G4double src[]) const;

  private:

    G4EquationOfMotion* fEquation = nullptr;
    G4int fnvar = 8;
    static constexpr G4double c_l = CLHEP::c_light/CLHEP::m*CLHEP::second;
};

#include "G4BorisScheme.icc"
#endif
