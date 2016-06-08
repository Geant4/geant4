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
// $Id: G4AlphaEvaporationChannel.cc,v 1.2.2.1 2001/06/28 19:13:10 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4AlphaEvaporationChannel.hh"


const G4AlphaEvaporationChannel & G4AlphaEvaporationChannel::operator=(const G4AlphaEvaporationChannel & right)
{
  G4Exception("G4AlphaEvaporationChannel::operator= meant to not be accessable");
  return *this;
}

G4AlphaEvaporationChannel::G4AlphaEvaporationChannel(const G4AlphaEvaporationChannel & right)
{
  G4Exception("G4AlphaEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4AlphaEvaporationChannel::operator==(const G4AlphaEvaporationChannel & right) const 
{
  return (this == (G4AlphaEvaporationChannel *) &right);
  //  return false;
}

G4bool G4AlphaEvaporationChannel::operator!=(const G4AlphaEvaporationChannel & right) const 
{
  return (this != (G4AlphaEvaporationChannel *) &right);
  //  return true;
}

