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
// $Id: G4DeuteronCoulombBarrier.cc,v 1.5 2008/09/19 13:32:54 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4DeuteronCoulombBarrier.hh"

G4DeuteronCoulombBarrier::G4DeuteronCoulombBarrier(const G4DeuteronCoulombBarrier & ) : G4CoulombBarrier()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4DeuteronCoulombBarrier & G4DeuteronCoulombBarrier::operator=(const G4DeuteronCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronCoulombBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4DeuteronCoulombBarrier::operator==(const G4DeuteronCoulombBarrier & ) const 
{
    return false;
}

G4bool G4DeuteronCoulombBarrier::operator!=(const G4DeuteronCoulombBarrier & ) const 
{
    return true;
}


G4double G4DeuteronCoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // const G4double Zlist[size] = {10.0, 20.0, 30.0, 50.0, 70.0};
    // const G4double Kprot[size] = {0.42, 0.58, 0.68, 0.77, 0.80};
    // 
    // K for deuteron is K for protons + 0.06
    G4double K = 1.0;
    if (aZ>=70.0) {
	K = 0.80;
    } else {
	K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
    }
    return K+0.06;
}
