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
// $Id: He3Channel.hh,v 1.1 2003-10-08 12:32:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef He3EvaporationChannel_h
#define He3EvaporationChannel_h 1

#include "EvapChannel.hh"
#include "G4He3CoulombBarrier.hh"
#include "G4He3EvaporationProbability.hh"

class He3EvaporationChannel : public EvaporationChannel
{
public:
  // only available constructor
  He3EvaporationChannel() : EvaporationChannel(3,2,"He3",
						   &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~He3EvaporationChannel() {};

private:
  const He3EvaporationChannel & operator=(const He3EvaporationChannel & right);  

  He3EvaporationChannel(const He3EvaporationChannel & right);

public:
  G4bool operator==(const He3EvaporationChannel & right) const;
  G4bool operator!=(const He3EvaporationChannel & right) const;

private:

  G4He3CoulombBarrier theCoulombBarrier;
	
  G4He3EvaporationProbability theEvaporationProbability;

};
#endif
