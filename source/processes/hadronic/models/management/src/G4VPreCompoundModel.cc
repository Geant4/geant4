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
// $Id: G4VPreCompoundModel.cc,v 1.3 2003/11/03 17:54:18 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#include "G4VPreCompoundModel.hh"
#include "G4HadronicException.hh"

G4VPreCompoundModel::G4VPreCompoundModel(G4ExcitationHandler *const value):
  theExcitationHandler(value)
{
}


const G4VPreCompoundModel & 
G4VPreCompoundModel::operator=(const G4VPreCompoundModel &)
{
  G4String text = "G4VPreCompoundModel::operator= meant to not be accessable";
  throw G4HadronicException(__FILE__, __LINE__, text); 
  return *this;
}

G4bool G4VPreCompoundModel::operator==(const G4VPreCompoundModel &) const
{
  return false;
}

G4bool G4VPreCompoundModel::operator!=(const G4VPreCompoundModel &) const
{
  return true;
}
