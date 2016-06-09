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
// $Id: G4TritonEvaporationChannel.cc,v 1.3 2005/06/04 13:21:21 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4TritonEvaporationChannel.hh"


const G4TritonEvaporationChannel & G4TritonEvaporationChannel::
operator=(const G4TritonEvaporationChannel & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonEvaporationChannel::operator= meant to not be accessable");
    return *this;
}

G4TritonEvaporationChannel::G4TritonEvaporationChannel(const G4TritonEvaporationChannel & ) : G4EvaporationChannel()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool G4TritonEvaporationChannel::operator==(const G4TritonEvaporationChannel & right) const 
{
    return (this == (G4TritonEvaporationChannel *) &right);
}

G4bool G4TritonEvaporationChannel::operator!=(const G4TritonEvaporationChannel & right) const 
{
    return (this != (G4TritonEvaporationChannel *) &right);
}
