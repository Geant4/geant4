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
// $Id: G4NeutronGEMChannel.hh,v 1.1 2003/08/26 18:43:04 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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
