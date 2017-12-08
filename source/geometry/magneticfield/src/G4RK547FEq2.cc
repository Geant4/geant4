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
//   0   |
//  2/13 |    2/13
//  2/13 |    3/52	          9/52
//  5/9  |    12955/26244     -15925/8748    12350/6561
//  3/4  |    -10383/52480    13923/10496    -176553/199424        505197/997120
//   1   |    1403/7236       -429/268       733330/309339         -7884/8911        104960/113967
//   1   |    181/2700        0              656903/1846800        19683/106400      34112/110565      67/800
//----------------------------------------------------------------------------------------------------------------------
//            181/2700        0              656903/1846800        19683/106400      34112/110565      67/800       0
//            11377/154575    0              35378291/105729300    343359/1522850    535952/1947645    134/17175    1/12

#include "G4RK547FEq2.hh"
#include "G4LineSection.hh"
#include "G4FieldUtils.hh"

using namespace field_utils;

namespace {

void copyArray(G4double dst[], const G4double src[])
{
    memcpy(dst, src, sizeof(G4double) * G4FieldTrack::ncompSVEC);
}

} // namespace

G4RK547FEq2::G4RK547FEq2(G4EquationOfMotion* EqRhs, G4int integrationVariables)
   : G4MagIntegratorStepper(EqRhs, integrationVariables)
{
}

void G4RK547FEq2::makeStep(
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
        b21 = 2./13.,
        b31 = 3./52., b32 = 9./52.,
        b41 = 12955./26244., b42 =  -15925./8748., b43 = 12350./6561.,
        b51 = -10383./52480., b52 = 13923./10496., b53 = -176553./199424.,
            b54 = 505197./997120.,
        b61 = 1403./7236., b62 = -429./268., b63 = 733330./309339.,
            b64 = -7884./8911., b65 = 104960./113967.,
        b71 = 181./2700., b72 = 0., b73 = 656903./1846800.,
             b74 = 19683./106400., b75 = 34112./110565., b76 =  67./800.;

    const G4double
        dc1 = b71 - 11377./154575.,
        dc2 = b72 - 0.,
        dc3 = b73 - 35378291./105729300.,
        dc4 = b74 - 343359./1522850.,
        dc5 = b75 - 535952./1947645.,
        dc6 = b76 - 134./17175.,
        dc7 = 0. - 1./12.;

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

void G4RK547FEq2::Stepper(
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

void G4RK547FEq2::Stepper(
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

G4double G4RK547FEq2::DistChord() const
{
    G4double yMid[G4FieldTrack::ncompSVEC];
    makeStep(fyIn, fdydx, fhstep / 2., yMid);

    const G4ThreeVector begin = makeVector(fyIn, Value3D::Position);
    const G4ThreeVector mid = makeVector(yMid, Value3D::Position);
    const G4ThreeVector end = makeVector(fyOut, Value3D::Position);

    return G4LineSection::Distline(mid, begin, end);
}
