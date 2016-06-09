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
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
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


