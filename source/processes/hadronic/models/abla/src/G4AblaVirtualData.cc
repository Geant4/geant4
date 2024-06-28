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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, UDC (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//

#include "G4AblaVirtualData.hh"
#include "globals.hh"

G4AblaVirtualData::G4AblaVirtualData() {}

G4bool G4AblaVirtualData::setAlpha(G4int A, G4int Z, G4double value)
{
    alpha[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setEcnz(G4int A, G4int Z, G4double value)
{
    ecnz[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setVgsld(G4int A, G4int Z, G4double value)
{
    vgsld[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setRms(G4int A, G4int Z, G4double value)
{
    rms[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setMexp(G4int A, G4int Z, G4double value)
{
    mexp[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setMexpID(G4int A, G4int Z, G4int value)
{
    mexpid[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setBeta2(G4int A, G4int Z, G4double value)
{
    beta2[A][Z] = value;

    return true;
}

G4bool G4AblaVirtualData::setBeta4(G4int A, G4int Z, G4double value)
{
    beta4[A][Z] = value;

    return true;
}

G4double G4AblaVirtualData::getAlpha(G4int A, G4int Z) { return alpha[A][Z]; }

G4double G4AblaVirtualData::getEcnz(G4int A, G4int Z) { return ecnz[A][Z]; }

G4double G4AblaVirtualData::getVgsld(G4int A, G4int Z) { return vgsld[A][Z]; }

G4double G4AblaVirtualData::getRms(G4int A, G4int Z) { return rms[A][Z]; }

G4double G4AblaVirtualData::getMexp(G4int A, G4int Z) { return mexp[A][Z]; }

G4int G4AblaVirtualData::getMexpID(G4int A, G4int Z) { return mexpid[A][Z]; }

G4double G4AblaVirtualData::getBeta2(G4int A, G4int Z) { return beta2[A][Z]; }

G4double G4AblaVirtualData::getBeta4(G4int A, G4int Z) { return beta4[A][Z]; }
