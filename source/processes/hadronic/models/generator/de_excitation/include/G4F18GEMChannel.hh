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
// $Id: G4F18GEMChannel.hh,v 1.1 2002/06/06 17:41:49 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4F18GEMChannel_h
#define G4F18GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4F18GEMCoulombBarrier.hh"
#include "G4F18GEMProbability.hh"

class G4F18GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4F18GEMChannel() : G4GEMChannel(18,9,"F18",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4F18GEMChannel() {};
  
private:
  const G4F18GEMChannel & operator=(const G4F18GEMChannel & right);  
    
  G4F18GEMChannel(const G4F18GEMChannel & right);
  
public:
  G4bool operator==(const G4F18GEMChannel & right) const;
  G4bool operator!=(const G4F18GEMChannel & right) const;
    
private:
  
  G4F18GEMCoulombBarrier theCoulombBarrier;
	
  G4F18GEMProbability theEvaporationProbability;
  
};
#endif
