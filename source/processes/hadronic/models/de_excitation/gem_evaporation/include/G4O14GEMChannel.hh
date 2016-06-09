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
// $Id: G4O14GEMChannel.hh,v 1.1 2003/08/26 18:43:05 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4O14GEMChannel_h
#define G4O14GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4O14GEMCoulombBarrier.hh"
#include "G4O14GEMProbability.hh"

class G4O14GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4O14GEMChannel() : G4GEMChannel(14,8,"O14",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4O14GEMChannel() {};
  
private:
  const G4O14GEMChannel & operator=(const G4O14GEMChannel & right);  
    
  G4O14GEMChannel(const G4O14GEMChannel & right);
  
public:
  G4bool operator==(const G4O14GEMChannel & right) const;
  G4bool operator!=(const G4O14GEMChannel & right) const;
    
private:
  
  G4O14GEMCoulombBarrier theCoulombBarrier;
	
  G4O14GEMProbability theEvaporationProbability;
  
};
#endif

