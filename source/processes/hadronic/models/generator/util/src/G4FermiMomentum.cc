// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FermiMomentum.cc,v 1.1 1998/08/22 08:58:07 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4FermiMomentum.hh"

G4FermiMomentum::G4FermiMomentum() : constofpmax(hbarc*cbrt(3.*pi2)) {}
G4FermiMomentum::~G4FermiMomentum(){}    
