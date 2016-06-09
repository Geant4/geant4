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
// $Id: G4DeuteronEvaporationChannel.hh,v 1.1 2003/08/26 18:30:41 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4DeuteronEvaporationChannel_h
#define G4DeuteronEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4DeuteronCoulombBarrier.hh"
#include "G4DeuteronEvaporationProbability.hh"

class G4DeuteronEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4DeuteronEvaporationChannel() : G4EvaporationChannel(2,1,"deuteron",
							&theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4DeuteronEvaporationChannel() {};

private:
  const G4DeuteronEvaporationChannel & operator=(const G4DeuteronEvaporationChannel & right);  

  G4DeuteronEvaporationChannel(const G4DeuteronEvaporationChannel & right);

public:
  G4bool operator==(const G4DeuteronEvaporationChannel & right) const;
  G4bool operator!=(const G4DeuteronEvaporationChannel & right) const;

private:

  G4DeuteronCoulombBarrier theCoulombBarrier;
	
  G4DeuteronEvaporationProbability theEvaporationProbability;

};
#endif
