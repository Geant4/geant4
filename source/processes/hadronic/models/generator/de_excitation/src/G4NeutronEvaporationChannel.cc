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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronEvaporationChannel.cc,v 1.4 2001/08/01 17:05:31 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4NeutronEvaporationChannel.hh"


const G4NeutronEvaporationChannel & G4NeutronEvaporationChannel::operator=(const G4NeutronEvaporationChannel & right)
{
    G4Exception("G4NeutronEvaporationChannel::operator= meant to not be accessable");
    return *this;
}

G4NeutronEvaporationChannel::G4NeutronEvaporationChannel(const G4NeutronEvaporationChannel & right)
{
    G4Exception("G4NeutronEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4NeutronEvaporationChannel::operator==(const G4NeutronEvaporationChannel & right) const 
{
    return (this == (G4NeutronEvaporationChannel *) &right);
    //  return false;
}

G4bool G4NeutronEvaporationChannel::operator!=(const G4NeutronEvaporationChannel & right) const 
{
    return (this != (G4NeutronEvaporationChannel *) &right);
    //  return true;
}

