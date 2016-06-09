//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ProtonCoulombBarrier.cc,v 1.2 2003/11/03 17:53:06 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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
