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
// $Id: G4Ne23GEMChannel.hh,v 1.1 2002/06/06 17:51:01 larazb Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4Ne23GEMChannel_h
#define G4Ne23GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4Ne23GEMCoulombBarrier.hh"
#include "G4Ne23GEMProbability.hh"

class G4Ne23GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4Ne23GEMChannel() : G4GEMChannel(23,10,"Ne23",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4Ne23GEMChannel() {};
  
private:
  const G4Ne23GEMChannel & operator=(const G4Ne23GEMChannel & right);  
    
  G4Ne23GEMChannel(const G4Ne23GEMChannel & right);
  
public:
  G4bool operator==(const G4Ne23GEMChannel & right) const;
  G4bool operator!=(const G4Ne23GEMChannel & right) const;
    
private:
  
  G4Ne23GEMCoulombBarrier theCoulombBarrier;
	
  G4Ne23GEMProbability theEvaporationProbability;
  
};
#endif
