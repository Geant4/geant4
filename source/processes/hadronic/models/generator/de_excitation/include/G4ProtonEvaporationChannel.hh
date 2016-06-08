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
// $Id: G4ProtonEvaporationChannel.hh,v 1.5 2002/12/12 19:17:11 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4ProtonEvaporationChannel_h
#define G4ProtonEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4ProtonCoulombBarrier.hh"
#include "G4ProtonEvaporationProbability.hh"

class G4ProtonEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4ProtonEvaporationChannel() : G4EvaporationChannel(1,1,"proton",
						      &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4ProtonEvaporationChannel() {};

private:
  const G4ProtonEvaporationChannel & operator=(const G4ProtonEvaporationChannel & right);  

  G4ProtonEvaporationChannel(const G4ProtonEvaporationChannel & right);

public:
  G4bool operator==(const G4ProtonEvaporationChannel & right) const;
  G4bool operator!=(const G4ProtonEvaporationChannel & right) const;


private:

  G4ProtonEvaporationProbability  theEvaporationProbability;
	
  G4ProtonCoulombBarrier  theCoulombBarrier;
};
#endif
