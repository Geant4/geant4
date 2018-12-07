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
//
// Description:
//    Modified midpoint method implementation
//    Implementation is based on modified_midpoint.hpp from boost odeint

//    Implementation by Dmitry Sorokin - GSoC 2016
//       Work supported by Google as part of Google Summer of Code 2016.
//    Supervision / code review: John Apostolakis
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4MODIFIED_MIDPOINT_HH
#define G4MODIFIED_MIDPOINT_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldTrack.hh"

class G4ModifiedMidpoint
{
  public:

    G4ModifiedMidpoint( G4EquationOfMotion* equation,
                        G4int nvar = 6, G4int steps = 2 );
   ~G4ModifiedMidpoint() = default;

    void DoStep( const G4double yIn[], const G4double dydxIn[],
                 G4double yOut[], G4double hstep) const;

    void DoStep( const G4double yIn[], const G4double dydxIn[],
                 G4double yOut[], G4double hstep, G4double yMid[],
                 G4double derivs[][G4FieldTrack::ncompSVEC]) const;

    inline void SetSteps(G4int steps);
    inline G4int GetSteps() const;

    inline void SetEquationOfMotion(G4EquationOfMotion* equation);
    inline G4EquationOfMotion* GetEquationOfMotion();

    inline G4int GetNumberOfVariables() const;

  private:

    void copy(G4double dst[], const G4double src[]) const;

  private:

    G4EquationOfMotion* fEquation;
    G4int fnvar;
    G4int fsteps;
};

#include "G4ModifiedMidpoint.icc"

#endif
