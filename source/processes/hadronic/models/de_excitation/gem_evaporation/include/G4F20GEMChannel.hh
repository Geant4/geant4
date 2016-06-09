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
// $Id: G4F20GEMChannel.hh,v 1.1 2003/08/26 18:42:08 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4F20GEMChannel_h
#define G4F20GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4F20GEMCoulombBarrier.hh"
#include "G4F20GEMProbability.hh"

class G4F20GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4F20GEMChannel() : G4GEMChannel(20,9,"F20",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4F20GEMChannel() {};
  
private:
  const G4F20GEMChannel & operator=(const G4F20GEMChannel & right);  
    
  G4F20GEMChannel(const G4F20GEMChannel & right);
  
public:
  G4bool operator==(const G4F20GEMChannel & right) const;
  G4bool operator!=(const G4F20GEMChannel & right) const;
    
private:
  
  G4F20GEMCoulombBarrier theCoulombBarrier;
	
  G4F20GEMProbability theEvaporationProbability;
  
};
#endif
