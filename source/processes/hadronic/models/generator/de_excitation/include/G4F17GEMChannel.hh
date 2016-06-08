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
// $Id: G4F17GEMChannel.hh,v 1.1 2002/06/06 17:41:46 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4F17GEMChannel_h
#define G4F17GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4F17GEMCoulombBarrier.hh"
#include "G4F17GEMProbability.hh"

class G4F17GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4F17GEMChannel() : G4GEMChannel(17,9,"F17",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4F17GEMChannel() {};
  
private:
  const G4F17GEMChannel & operator=(const G4F17GEMChannel & right);  
    
  G4F17GEMChannel(const G4F17GEMChannel & right);
  
public:
  G4bool operator==(const G4F17GEMChannel & right) const;
  G4bool operator!=(const G4F17GEMChannel & right) const;
    
private:
  
  G4F17GEMCoulombBarrier theCoulombBarrier;
	
  G4F17GEMProbability theEvaporationProbability;
  
};
#endif
