// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4VPhotonEvaporation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 28 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------


#include "G4VPhotonEvaporation.hh"


G4bool G4VPhotonEvaporation::operator==(const G4VPhotonEvaporation &right) const
{
    return (this == (G4VPhotonEvaporation*) &right);
}

G4bool G4VPhotonEvaporation::operator!=(const G4VPhotonEvaporation &right) const
{
    return (this != (G4VPhotonEvaporation*) &right);
}


