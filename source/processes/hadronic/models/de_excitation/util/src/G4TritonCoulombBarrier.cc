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
// $Id: G4TritonCoulombBarrier.cc,v 1.3 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4TritonCoulombBarrier.hh"

G4TritonCoulombBarrier::G4TritonCoulombBarrier(const G4TritonCoulombBarrier & ) : G4CoulombBarrier()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4TritonCoulombBarrier & G4TritonCoulombBarrier::operator=(const G4TritonCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonCoulombBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4TritonCoulombBarrier::operator==(const G4TritonCoulombBarrier & ) const 
{
    return false;
}

G4bool G4TritonCoulombBarrier::operator!=(const G4TritonCoulombBarrier & ) const 
{
    return true;
}


G4double G4TritonCoulombBarrier::BarrierPenetrationFactor(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // const G4double Zlist[size] = {10.0, 20.0, 30.0, 50.0, 70.0};
    // const G4double Kprot[size] = {0.42, 0.58, 0.68, 0.77, 0.80};
    // 
    // K for Triton is K for protons + 0.12
    G4double K = 1.0;
    if (aZ>=70.0) {
	K = 0.80;
    } else {
	K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
    }
    return K+0.12;
}
