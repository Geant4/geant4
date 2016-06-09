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
// $Id: G4AlphaEvaporationChannel.hh,v 1.2 2005/06/04 13:21:21 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4AlphaEvaporationChannel_h
#define G4AlphaEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4AlphaCoulombBarrier.hh"
#include "G4AlphaEvaporationProbability.hh"

class G4AlphaEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4AlphaEvaporationChannel() : G4EvaporationChannel(4,2,"alpha",
						     &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4AlphaEvaporationChannel() {};

private:
  const G4AlphaEvaporationChannel & operator=(const G4AlphaEvaporationChannel & right);  

  G4AlphaEvaporationChannel(const G4AlphaEvaporationChannel & right);


public:
  G4bool operator==(const G4AlphaEvaporationChannel & right) const;
  G4bool operator!=(const G4AlphaEvaporationChannel & right) const;

private:

  G4AlphaCoulombBarrier theCoulombBarrier;
	
  G4AlphaEvaporationProbability theEvaporationProbability;

};
#endif
