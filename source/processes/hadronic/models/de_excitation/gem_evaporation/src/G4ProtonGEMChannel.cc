//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4ProtonGEMChannel.cc,v 1.5 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//

#include "G4ProtonGEMChannel.hh"


const G4ProtonGEMChannel & G4ProtonGEMChannel::operator=(const G4ProtonGEMChannel & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonGEMChannel::operator= meant to not be accessable");
    return *this;
}

G4ProtonGEMChannel::G4ProtonGEMChannel(const G4ProtonGEMChannel & ): G4GEMChannel()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonGEMChannel::CopyConstructor meant to not be accessable");
}

G4bool G4ProtonGEMChannel::operator==(const G4ProtonGEMChannel & right) const 
{
    return (this == (G4ProtonGEMChannel *) &right);
    //  return false;
}

G4bool G4ProtonGEMChannel::operator!=(const G4ProtonGEMChannel & right) const 
{
    return (this != (G4ProtonGEMChannel *) &right);
    //  return true;
}

