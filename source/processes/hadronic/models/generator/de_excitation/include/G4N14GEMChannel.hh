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
// $Id: G4N14GEMChannel.hh,v 1.1 2002/06/06 17:49:52 larazb Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4N14GEMChannel_h
#define G4N14GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4N14GEMCoulombBarrier.hh"
#include "G4N14GEMProbability.hh"

class G4N14GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4N14GEMChannel() : G4GEMChannel(14,7,"N14",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4N14GEMChannel() {};
  
private:
  const G4N14GEMChannel & operator=(const G4N14GEMChannel & right);  
    
  G4N14GEMChannel(const G4N14GEMChannel & right);
  
public:
  G4bool operator==(const G4N14GEMChannel & right) const;
  G4bool operator!=(const G4N14GEMChannel & right) const;
    
private:
  
  G4N14GEMCoulombBarrier theCoulombBarrier;
	
  G4N14GEMProbability theEvaporationProbability;
  
};
#endif
