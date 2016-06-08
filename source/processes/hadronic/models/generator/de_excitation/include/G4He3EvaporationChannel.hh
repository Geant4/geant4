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
// $Id: G4He3EvaporationChannel.hh,v 1.5 2002/12/12 19:17:06 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4He3EvaporationChannel_h
#define G4He3EvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4He3CoulombBarrier.hh"
#include "G4He3EvaporationProbability.hh"

class G4He3EvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4He3EvaporationChannel() : G4EvaporationChannel(3,2,"He3",
						   &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4He3EvaporationChannel() {};

private:
  const G4He3EvaporationChannel & operator=(const G4He3EvaporationChannel & right);  

  G4He3EvaporationChannel(const G4He3EvaporationChannel & right);

public:
  G4bool operator==(const G4He3EvaporationChannel & right) const;
  G4bool operator!=(const G4He3EvaporationChannel & right) const;

private:

  G4He3CoulombBarrier theCoulombBarrier;
	
  G4He3EvaporationProbability theEvaporationProbability;

};
#endif
