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
// $Id: G4NeutronGEMChannel.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4NeutronGEMChannel_h
#define G4NeutronGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4NeutronCoulombBarrier.hh"
#include "G4NeutronGEMProbability.hh"

class G4NeutronGEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4NeutronGEMChannel() : G4GEMChannel(1,0,"neutron",
                                       &theEvaporationProbability,
                                       &theCoulombBarrier) {};
  // destructor
  ~G4NeutronGEMChannel() {};

private:
  const G4NeutronGEMChannel & operator=(const G4NeutronGEMChannel & right);  

  G4NeutronGEMChannel(const G4NeutronGEMChannel & right);

public:
  G4bool operator==(const G4NeutronGEMChannel & right) const;
  G4bool operator!=(const G4NeutronGEMChannel & right) const;

private:

  G4NeutronCoulombBarrier theCoulombBarrier;
	
  G4NeutronGEMProbability theEvaporationProbability;

};
#endif
