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
// $Id: G4O19GEMChannel.hh,v 1.1 2002/06/06 17:52:13 larazb Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4O19GEMChannel_h
#define G4O19GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4O19GEMCoulombBarrier.hh"
#include "G4O19GEMProbability.hh"

class G4O19GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4O19GEMChannel() : G4GEMChannel(19,8,"O19",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4O19GEMChannel() {};
  
private:
  const G4O19GEMChannel & operator=(const G4O19GEMChannel & right);  
    
  G4O19GEMChannel(const G4O19GEMChannel & right);
  
public:
  G4bool operator==(const G4O19GEMChannel & right) const;
  G4bool operator!=(const G4O19GEMChannel & right) const;
    
private:
  
  G4O19GEMCoulombBarrier theCoulombBarrier;
	
  G4O19GEMProbability theEvaporationProbability;
  
};
#endif

