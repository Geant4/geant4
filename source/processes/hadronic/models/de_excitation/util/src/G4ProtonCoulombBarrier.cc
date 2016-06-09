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
// $Id: G4ProtonCoulombBarrier.cc,v 1.4 2006/06/29 20:28:33 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4ProtonCoulombBarrier.hh"

G4ProtonCoulombBarrier::G4ProtonCoulombBarrier(const G4ProtonCoulombBarrier & ) : G4CoulombBarrier()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4ProtonCoulombBarrier & G4ProtonCoulombBarrier::operator=(const G4ProtonCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonCoulombBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4ProtonCoulombBarrier::operator==(const G4ProtonCoulombBarrier & ) const 
{
    return false;
}

G4bool G4ProtonCoulombBarrier::operator!=(const G4ProtonCoulombBarrier & ) const 
{
    return true;
}


G4double G4ProtonCoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // const G4double Zlist[size] = {10.0, 20.0, 30.0, 50.0, 70.0};
    // const G4double Kprot[size] = {0.42, 0.58, 0.68, 0.77, 0.80};
    G4double K = 1.0;
    if (aZ>=70.0) {
	K = 0.80;
    } else {
	K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
    }
    return K;
}
