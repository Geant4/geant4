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
// $Id: AlphaChannel.cc,v 1.1 2003-10-08 12:32:18 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "AlphaChannel.hh"


const AlphaEvaporationChannel & AlphaEvaporationChannel::operator=(const AlphaEvaporationChannel & right)
{
  G4Exception("G4AlphaEvaporationChannel::operator= meant to not be accessable");
  return *this;
}

AlphaEvaporationChannel::AlphaEvaporationChannel(const AlphaEvaporationChannel & right)
{
  G4Exception("G4AlphaEvaporationChannel::CopyConstructor meant to not be accessable");
}

G4bool AlphaEvaporationChannel::operator==(const AlphaEvaporationChannel & right) const 
{
  return (this == (AlphaEvaporationChannel *) &right);
  //  return false;
}

G4bool AlphaEvaporationChannel::operator!=(const AlphaEvaporationChannel & right) const 
{
  return (this != (AlphaEvaporationChannel *) &right);
  //  return true;
}

