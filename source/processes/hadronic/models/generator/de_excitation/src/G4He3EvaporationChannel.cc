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
// $Id: G4He3EvaporationChannel.cc,v 1.2.2.1 2001/06/28 19:13:17 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4He3EvaporationChannel.hh"


const G4He3EvaporationChannel & G4He3EvaporationChannel::operator=(const G4He3EvaporationChannel & right)
{
    G4Exception("G4He3EvaporationChannel::operator= meant to not be accessable");
    return *this;
}

G4He3EvaporationChannel::G4He3EvaporationChannel(const G4He3EvaporationChannel & right)
{
    G4Exception("G4He3EvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4He3EvaporationChannel::operator==(const G4He3EvaporationChannel & right) const 
{
    return (this == (G4He3EvaporationChannel *) &right);
    //  return false;
}

G4bool G4He3EvaporationChannel::operator!=(const G4He3EvaporationChannel & right) const 
{
    return (this != (G4He3EvaporationChannel *) &right);
    //  return true;
}

