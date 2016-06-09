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
// $Id: G4TritonGEMChannel.cc,v 1.3 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4TritonGEMChannel.hh"


const G4TritonGEMChannel & G4TritonGEMChannel::operator=(const G4TritonGEMChannel & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonGEMChannel::operator= meant to not be accessable");
    return *this;
}

G4TritonGEMChannel::G4TritonGEMChannel(const G4TritonGEMChannel & ) : G4GEMChannel()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonGEMChannel::CopyConstructor meant to not be accessable");
}

G4bool G4TritonGEMChannel::operator==(const G4TritonGEMChannel & right) const 
{
    return (this == (G4TritonGEMChannel *) &right);
    //  return false;
}

G4bool G4TritonGEMChannel::operator!=(const G4TritonGEMChannel & right) const 
{
    return (this != (G4TritonGEMChannel *) &right);
    //  return true;
}

