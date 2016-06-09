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
// $Id: G4O15GEMChannel.hh,v 1.1 2003/08/26 18:43:07 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4O15GEMChannel_h
#define G4O15GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4O15GEMCoulombBarrier.hh"
#include "G4O15GEMProbability.hh"

class G4O15GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4O15GEMChannel() : G4GEMChannel(15,8,"O15",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4O15GEMChannel() {};
  
private:
  const G4O15GEMChannel & operator=(const G4O15GEMChannel & right);  
    
  G4O15GEMChannel(const G4O15GEMChannel & right);
  
public:
  G4bool operator==(const G4O15GEMChannel & right) const;
  G4bool operator!=(const G4O15GEMChannel & right) const;
    
private:
  
  G4O15GEMCoulombBarrier theCoulombBarrier;
	
  G4O15GEMProbability theEvaporationProbability;
  
};
#endif

