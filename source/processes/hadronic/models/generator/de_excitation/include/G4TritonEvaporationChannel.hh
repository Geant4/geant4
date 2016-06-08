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
// $Id: G4TritonEvaporationChannel.hh,v 1.4 2001/08/01 17:04:52 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4TritonEvaporationChannel_h
#define G4TritonEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4TritonCoulombBarrier.hh"
#include "G4TritonEvaporationProbability.hh"

class G4TritonEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4TritonEvaporationChannel() : G4EvaporationChannel(3,1,"triton",
						      &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4TritonEvaporationChannel() {};

private:
  const G4TritonEvaporationChannel & operator=(const G4TritonEvaporationChannel & right);  

  G4TritonEvaporationChannel(const G4TritonEvaporationChannel & right);

public:
  G4bool operator==(const G4TritonEvaporationChannel & right) const;
  G4bool operator!=(const G4TritonEvaporationChannel & right) const;

private:

  G4TritonCoulombBarrier theCoulombBarrier;
	
  G4TritonEvaporationProbability theEvaporationProbability;

};
#endif
