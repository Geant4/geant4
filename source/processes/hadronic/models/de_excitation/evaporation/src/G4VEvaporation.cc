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
// $Id: G4VEvaporation.cc,v 1.4 2005/06/04 13:21:21 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) writen from G4Evaporation.cc (May 1998)
//


#include "G4VEvaporation.hh"
#include "G4HadronicException.hh"

G4VEvaporation::G4VEvaporation(const G4VEvaporation &)
{
 throw G4HadronicException(__FILE__, __LINE__, "G4VEvaporation::copy_constructor meant to not be accessable");
}




const G4VEvaporation & G4VEvaporation::operator=(const G4VEvaporation &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4VEvaporation::operator= meant to not be accessable");
  return *this;
}


G4bool G4VEvaporation::operator==(const G4VEvaporation &) const
{
  return false;
}

G4bool G4VEvaporation::operator!=(const G4VEvaporation &) const
{
  return true;
}





