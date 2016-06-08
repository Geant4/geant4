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
// $Id: G4NeutronEvaporationChannel.hh,v 1.4 2001/08/01 17:04:34 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4NeutronEvaporationChannel_h
#define G4NeutronEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4NeutronCoulombBarrier.hh"
#include "G4NeutronEvaporationProbability.hh"

class G4NeutronEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4NeutronEvaporationChannel() : G4EvaporationChannel(1,0,"neutron",
						       &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4NeutronEvaporationChannel() {};

private:
  const G4NeutronEvaporationChannel & operator=(const G4NeutronEvaporationChannel & right);  

  G4NeutronEvaporationChannel(const G4NeutronEvaporationChannel & right);

public:
  G4bool operator==(const G4NeutronEvaporationChannel & right) const;
  G4bool operator!=(const G4NeutronEvaporationChannel & right) const;

private:

  G4NeutronCoulombBarrier theCoulombBarrier;
	
  G4NeutronEvaporationProbability theEvaporationProbability;

};
#endif
