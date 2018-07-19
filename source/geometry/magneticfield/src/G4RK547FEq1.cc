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
// The Butcher table of the Higham & Hall 5(4)7 method is:
//
//   0  |
//  2/9 |      2/9
//  1/3 |      1/12	     1/4
//  1/2 |      1/8       0          3/8
//  3/5 |      91/500    -27/100    78/125    8/125
//   1  |      -11/20    27/20      12/5      -36/5    5
//   1  |      1/12      0          27/32     -4/3     125/96    5/48
//----------------------------------------------------------------------------
//             1/12      0          27/32     -4/3     125/96    5/48    0
//             2/15      0          27/80     -2/15    25/48     1/24    1/10

#include "G4RK547FEq1.hh"
#include "G4LineSection.hh"
#include "G4FieldUtils.hh"

using namespace field_utils;

namespace {

void copyArray(G4double dst[], const G4double src[])
{
    memcpy(dst, src, sizeof(G4double) * G4FieldTrack::ncompSVEC);
}

} // namespace

G4RK547FEq1::G4RK547FEq1(G4EquationOfMotion* EqRhs, G4int integrationVariables)
   : G4MagIntegratorStepper(EqRhs, integrationVariables)
{
}

void G4RK547FEq1::makeStep(
    const G4double yInput[],
    const G4double dydx[],
    const G4double hstep,
    G4double yOutput[],
    G4double* dydxOutput,
    G4double* yError) const
{
    G4double yTemp[G4FieldTrack::ncompSVEC];
    for (int i = GetNumberOfVariables(); i < GetNumberOfStateVariables(); ++i){
        yOutput[i] = yTemp[i] = yInput[i];
    }

    G4double ak2[G4FieldTrack::ncompSVEC],
             ak3[G4FieldTrack::ncompSVEC],
             ak4[G4FieldTrack::ncompSVEC],
             ak5[G4FieldTrack::ncompSVEC],
             ak6[G4FieldTrack::ncompSVEC];

    const G4double
        b21 = 2./9.,
        b31 = 1./12., b32 = 1./4.,
        b41 = 1./8., b42 = 0., b43 = 3./8.,
        b51 = 91./500., b52 = -27./100., b53 = 78./125., b54 = 8./125.,
        b61 = -11./20., b62 = 27./20., b63 = 12./5.,
            b64 = -36./5., b65 = 5.,
        b71 = 1./12.,    b72 = 0., b73 = 27./32.,
             b74 = -4./3., b75 = 125./96., b76 = 5./48.;

    const G4double
        dc1 = b71 - 2./15.,
        dc2 = b72 - 0.,
        dc3 = b73 - 27./80.,
        dc4 = b74 + 2./15.,
        dc5 = b75 - 25./48.,
        dc6 = b76 - 1./24.,
        dc7 = 0. - 1./10.;

    //RightHandSide(yInput, dydx);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp[i] = yInput[i] + hstep * b21 * dydx[i];

    RightHandSide(yTemp, ak2);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp[i] = yInput[i] + hstep * (b31 * dydx[i] + b32 * ak2[i]);

    RightHandSide(yTemp, ak3);
    for(int i = 0;i < GetNumberOfVariables(); ++i)
        yTemp[i] = yInput[i] + hstep * (b41 * dydx[i] + b42 * ak2[i] +
                                        b43 * ak3[i]);

    RightHandSide(yTemp, ak4);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp[i] = yInput[i] + hstep * (b51 * dydx[i] + b52 * ak2[i] +
                                        b53 * ak3[i] + b54 * ak4[i]);

    RightHandSide(yTemp, ak5);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp[i] = yInput[i] + hstep * (b61 * dydx[i] + b62 * ak2[i] +
                                        b63 * ak3[i] + b64 * ak4[i] +
                                        b65 * ak5[i]);

    RightHandSide(yTemp, ak6);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yOutput[i] = yInput[i] + hstep * (b71 * dydx[i] + b72 * ak2[i] +
                                          b73 * ak3[i] + b74 * ak4[i] +
                                          b75 * ak5[i] + b76 * ak6[i]);

    if (dydxOutput && yError) {
        RightHandSide(yOutput, dydxOutput);
        for(int i = 0; i < GetNumberOfVariables(); ++i)
            yError[i] = hstep * (dc1 * dydx[i] + dc2 * ak2[i] + dc3 * ak3[i] +
                                 dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i] +
                                 dc7 * dydxOutput[i]);
    }
}

void G4RK547FEq1::Stepper(
    const G4double yInput[],
    const G4double dydx[],
    G4double hstep,
    G4double yOutput[],
    G4double yError[])
{
    copyArray(fyIn, yInput);
    copyArray(fdydx, dydx);
    fhstep = hstep;

    makeStep(fyIn, fdydx, fhstep, fyOut, fdydxOut, yError);

    copyArray(yOutput, fyOut);
}

void G4RK547FEq1::Stepper(
    const G4double yInput[],
    const G4double dydx[],
    G4double hstep,
    G4double yOutput[],
    G4double yError[],
    G4double dydxOutput[])
{
    copyArray(fyIn, yInput);
    copyArray(fdydx, dydx);
    fhstep = hstep;

    makeStep(fyIn,fdydx, fhstep, fyOut, fdydxOut, yError);

    copyArray(yOutput, fyOut);
    copyArray(dydxOutput, fdydxOut);
}

G4double G4RK547FEq1::DistChord() const
{
    G4double yMid[G4FieldTrack::ncompSVEC];
    makeStep(fyIn, fdydx, fhstep / 2., yMid);

    const G4ThreeVector begin = makeVector(fyIn, Value3D::Position);
    const G4ThreeVector mid = makeVector(yMid, Value3D::Position);
    const G4ThreeVector end = makeVector(fyOut, Value3D::Position);

    return G4LineSection::Distline(mid, begin, end);
}
