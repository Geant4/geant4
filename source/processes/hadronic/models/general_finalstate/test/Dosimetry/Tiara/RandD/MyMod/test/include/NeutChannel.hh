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
// $Id: NeutChannel.hh,v 1.1 2003-10-08 12:32:18 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef NeutronEvaporationChannel_h
#define NeutronEvaporationChannel_h 1

#include "EvapChannel.hh"
#include "G4NeutronCoulombBarrier.hh"
#include "G4NeutronEvaporationProbability.hh"

class NeutronEvaporationChannel : public EvaporationChannel
{
public:
  // only available constructor
  NeutronEvaporationChannel() : EvaporationChannel(1,0,"neutron",
						       &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~NeutronEvaporationChannel() {};

private:
  const NeutronEvaporationChannel & operator=(const NeutronEvaporationChannel & right);  

  NeutronEvaporationChannel(const NeutronEvaporationChannel & right);

public:
  G4bool operator==(const NeutronEvaporationChannel & right) const;
  G4bool operator!=(const NeutronEvaporationChannel & right) const;

private:

  G4NeutronCoulombBarrier theCoulombBarrier;
	
  G4NeutronEvaporationProbability theEvaporationProbability;

};
#endif
