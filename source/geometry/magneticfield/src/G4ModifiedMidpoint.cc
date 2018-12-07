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
//    G4ModifiedMidpoint implementation
//    Based on modified_midpoint.hpp from boost
//
//    Author: Dmitry Sorokin - GSoC 2016
//
///////////////////////////////////////////////////////////////////////////////

#include "G4ModifiedMidpoint.hh"
#include "G4FieldUtils.hh"

using namespace field_utils;

G4ModifiedMidpoint::G4ModifiedMidpoint( G4EquationOfMotion* equation,
                                        G4int nvar, G4int steps )
  : fEquation(equation), fnvar(nvar), fsteps(steps)
{
  if (nvar <= 0)
  {
    G4Exception("G4ModifiedMidpoint::G4ModifiedMidpoint()",
                "GeomField0002", FatalException,
                "Invalid number of variables; must be greater than zero!");
  }
}

void G4ModifiedMidpoint::DoStep( const G4double yIn[], const G4double dydyIn[],
                                 G4double yOut[], G4double hstep) const
{
  G4double y0[G4FieldTrack::ncompSVEC];
  G4double y1[G4FieldTrack::ncompSVEC];
  G4double yTemp[G4FieldTrack::ncompSVEC];
  setValue(yIn, Value1D::LabTime, y0, y1, yTemp, yOut);

  G4double dydx[G4FieldTrack::ncompSVEC];

  const G4double h = hstep / fsteps;
  const G4double h2 = 2 * h;

  // y1 = yIn + h * dydx
  for (G4int i = 0; i < fnvar; ++i)
  {
    y1[i] = yIn[i] + h * dydyIn[i];
  }

  fEquation->RightHandSide(y1, dydx);

  copy(y0, yIn);

  // general step
  // yTemp = y1; y1 = y0 + h2 * dydx; y0 = yTemp
  for (G4int i = 1; i < fsteps; ++i)
  {
    copy(yTemp, y1);
    for (G4int j = 0; j < fnvar; ++j)
    {
      y1[j] = y0[j] + h2 * dydx[j];
    }
    copy(y0, yTemp);

    fEquation->RightHandSide(y1, dydx);
  }

  // last step
  // yOut = 0.5 * (y0 + y1 + h * dydx)
  for (G4int i = 0; i < fnvar; ++i)
  {
    yOut[i] = 0.5 * (y0[i] + y1[i] + h * dydx[i]);
  }
}

void G4ModifiedMidpoint::DoStep( const G4double yIn[], const G4double dydxIn[],
                           G4double yOut[], G4double hstep, G4double yMid[],
                           G4double derivs[][G4FieldTrack::ncompSVEC]) const
{
  G4double y0[G4FieldTrack::ncompSVEC];
  G4double y1[G4FieldTrack::ncompSVEC];
  G4double yTemp[G4FieldTrack::ncompSVEC];
  setValue(yIn, Value1D::LabTime, y0, y1, yTemp, yMid, yOut);

  const G4double h = hstep / fsteps;
  const G4double h2 = 2 * h;

  // y0 = yIn
  copy(y0, yIn);

  // y1 = y0 + h * dydx
  for (G4int i = 0; i < fnvar; ++i)
  {
    y1[i] = y0[i] + h * dydxIn[i];
  }

  // result of first step already gives approximation
  // at the center of the interval
  if(fsteps == 2)
  {
    copy(yMid, y1);
  }

  fEquation->RightHandSide(y1, derivs[0]);

  // general step
  // yTemp = y1; y1 = y0 + h2 * dydx; y0 = yTemp
  for (G4int i = 1; i < fsteps; ++i)
  {
    copy(yTemp, y1);
    for (G4int j = 0; j < fnvar; ++j)
    {
      y1[j] = y0[j] + h2 * derivs[i-1][j];
    }
    copy(y0, yTemp);

    // save approximation at the center of the interval
    if(i == fsteps / 2 - 1 )
    {
      copy(yMid, y1);
    }

    fEquation->RightHandSide(y1, derivs[i]);
  }

  // last step
  // yOut = 0.5 * (y0 + y1 + h * dydx)
  for (G4int i = 0; i < fnvar; ++i)
  {
    yOut[i] = 0.5 * (y0[i] + y1[i] + h * derivs[fsteps-1][i]);
  }
}

void G4ModifiedMidpoint::copy(G4double dst[], const G4double src[]) const
{
  std::memcpy(dst, src, sizeof(G4double) * fnvar);
}
