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
//   0    |
//  11/45 |    11/45
//  11/30 |    11/120          11/40
//  55/56 |    106865/87808    -408375/87808    193875/43904
//  9/10  |    79503/121000    -1053/440        147753/56870       27048/710875
//   1    |    89303/78045     -2025/473        994650/244541      -2547216/28122215     475/2967
//   1    |    1247/10890      0                57375/108053       -1229312/1962015      125/207     43/114
//---------------------------------------------------------------------------------------------------------------------
//             1247/10890      0                57375/108053       -1229312/1962015      125/207     43/114       0
//             21487/185130    0                963225/1836901     -39864832/33354255    2575/3519   4472/4845    -1/10

#include "G4RK547FEq3.hh"
#include "G4LineSection.hh"
#include "G4FieldUtils.hh"

using namespace field_utils;

namespace {

void copyArray(G4double dst[], const G4double src[])
{
    memcpy(dst, src, sizeof(G4double) * G4FieldTrack::ncompSVEC);
}

} // namespace

G4RK547FEq3::G4RK547FEq3(G4EquationOfMotion* EqRhs, G4int integrationVariables)
   : G4MagIntegratorStepper(EqRhs, integrationVariables)
{
}

void G4RK547FEq3::makeStep(
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
        b21 = 11./45.,
        b31 = 11./120., b32 = 11./40.,
        b41 = 106865./87808., b42 = -408375./87808., b43 = 193875./43904.,
        b51 = 79503./121000., b52 = -1053./440., b53 = 147753./56870.,
            b54 =  27048./710875.,
        b61 = 89303./78045., b62 = -2025./473., b63 = 994650./244541.,
            b64 = -2547216./28122215., b65 = 475./2967.,
        b71 = 1247./10890., b72 = 0., b73 = 57375./108053.,
             b74 = -1229312./1962015., b75 = 125./207., b76 = 43./114.;

    const G4double
        dc1 = b71 - 21487./185130.,
        dc2 = b72 - 0.,
        dc3 = b73 - 963225./1836901.,
        dc4 = b74 + 39864832./33354255.,
        dc5 = b75 - 2575./3519.,
        dc6 = b76 - 4472./4845.,
        dc7 = 0. + 1./10.;

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

void G4RK547FEq3::Stepper(
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

void G4RK547FEq3::Stepper(
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

    makeStep(fyIn, fdydx, fhstep, fyOut, fdydxOut, yError);

    copyArray(yOutput, fyOut);
    copyArray(dydxOutput, fdydxOut);
}

G4double G4RK547FEq3::DistChord() const
{
    G4double yMid[G4FieldTrack::ncompSVEC];
    makeStep(fyIn, fdydx, fhstep / 2., yMid);

    const G4ThreeVector begin = makeVector(fyIn, Value3D::Position);
    const G4ThreeVector mid = makeVector(yMid, Value3D::Position);
    const G4ThreeVector end = makeVector(fyOut, Value3D::Position);

    return G4LineSection::Distline(mid, begin, end);
}
