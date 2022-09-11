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
// G4BorisSDC
//
// Class description:
//
// Implimentation of the boris algorithm for advancing 
// charged particles in an elctromagnetic field.

// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun
// --------------------------------------------------------------------
#ifndef G4BORIS_SDC_HH
#define G4BORIS_SDC_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldTrack.hh"
#include "G4PhysicalConstants.hh"
#include"G4SystemOfUnits.hh"

class G4BorisSDC
{
  public:

    G4BorisSDC() = default;
    G4BorisSDC( G4EquationOfMotion* equation,
                        G4int nvar = 6);
   ~G4BorisSDC() = default;

    G4ThreeVector GetLorentzForce(const G4ThreeVector position, const G4ThreeVector velocity) const;

    void DoStep(const G4double restMass, const G4double charge, const G4double yIn[], 
                 G4double yOut[], G4double hstep);

    void UpdatePosition(G4int k, G4int m) ;

    void UpdateVelocity( G4int k, G4int m) ;
    
    

    

    inline void SetEquationOfMotion(G4EquationOfMotion* equation);
    inline G4EquationOfMotion* GetEquationOfMotion();

    inline G4int GetNumberOfVariables() const;

  private:

    void copy(G4double dst[], const G4double src[]) const;

  private:

    G4EquationOfMotion* fEquation = nullptr;
    G4int fnvar = 0;
    static constexpr G4double c_l = CLHEP::c_light/CLHEP::m*CLHEP::second;
    G4double alpha; // charge/mass ratio (SI)
    G4double mass_si;
    G4double restMass_c2;
    G4double charge_si;
    static constexpr G4int M = 3; // no. of nodes(3) for the order 2M-2 = 4
    static constexpr G4int K = 4; // no. of iterations 
    static constexpr double sqrt15 = 3.872983346207417;
    static constexpr G4double nodes[M + 1] ={0, 0.5 - sqrt15/10.0, 0.5, 0.5 + sqrt15/10.0 };
    G4double delta_t_m[M + 1];
    static constexpr G4double Q[M+1][M+1] ={0, 0, 0, 0, 
                           0, 5.0/36, 2.0/9 - sqrt15/15, 5.0/36 - sqrt15/30,
                           0, 5.0/36 + sqrt15/24, 2.0/9, 5.0/36 - sqrt15/24,
                           0, 5.0/36 + sqrt15/30, 2.0/9 + sqrt15/15, 5.0/36 };
    G4double S[M+1][M+1];
    // Initialize the starting values of velocity and position at various nodes and iterations
    G4ThreeVector Velocity[M+1][K+1]; 
    G4ThreeVector Position[M+1][K+1]; 

    
    

    
    
    
};

#include "G4BorisSDC.icc"
#endif