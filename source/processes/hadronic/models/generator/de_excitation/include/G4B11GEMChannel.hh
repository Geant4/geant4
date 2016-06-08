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
// $Id: G4B11GEMChannel.hh,v 1.1 2002/06/06 17:37:12 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4B11GEMChannel_h
#define G4B11GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4B11GEMCoulombBarrier.hh"
#include "G4B11GEMProbability.hh"

class G4B11GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4B11GEMChannel() : G4GEMChannel(11,5,"B11",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4B11GEMChannel() {};
  
private:
  const G4B11GEMChannel & operator=(const G4B11GEMChannel & right);  
    
  G4B11GEMChannel(const G4B11GEMChannel & right);
  
public:
  G4bool operator==(const G4B11GEMChannel & right) const;
  G4bool operator!=(const G4B11GEMChannel & right) const;
    
private:
  
  G4B11GEMCoulombBarrier theCoulombBarrier;
	
  G4B11GEMProbability theEvaporationProbability;
  
};
#endif
