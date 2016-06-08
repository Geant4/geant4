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
// $Id: G4Be10GEMChannel.hh,v 1.1 2002/06/06 17:39:25 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4Be10GEMChannel_h
#define G4Be10GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4Be10GEMCoulombBarrier.hh"
#include "G4Be10GEMProbability.hh"

class G4Be10GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4Be10GEMChannel() : G4GEMChannel(10,4,"Be10",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4Be10GEMChannel() {};
  
private:
  const G4Be10GEMChannel & operator=(const G4Be10GEMChannel & right);  
    
  G4Be10GEMChannel(const G4Be10GEMChannel & right);
  
public:
  G4bool operator==(const G4Be10GEMChannel & right) const;
  G4bool operator!=(const G4Be10GEMChannel & right) const;
    
private:
  
  G4Be10GEMCoulombBarrier theCoulombBarrier;
	
  G4Be10GEMProbability theEvaporationProbability;
  
};
#endif
