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
// $Id: G4He8GEMChannel.hh,v 1.1 2002/06/06 17:43:14 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4He8GEMChannel_h
#define G4He8GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4He8GEMCoulombBarrier.hh"
#include "G4He8GEMProbability.hh"

class G4He8GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4He8GEMChannel() : G4GEMChannel(8,2,"He8",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4He8GEMChannel() {};
  
private:
  const G4He8GEMChannel & operator=(const G4He8GEMChannel & right);  
    
  G4He8GEMChannel(const G4He8GEMChannel & right);
  
public:
  G4bool operator==(const G4He8GEMChannel & right) const;
  G4bool operator!=(const G4He8GEMChannel & right) const;
    
private:
  
  G4He8GEMCoulombBarrier theCoulombBarrier;
	
  G4He8GEMProbability theEvaporationProbability;
  
};
#endif
