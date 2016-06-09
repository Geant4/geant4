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
// $Id: G4Be9GEMChannel.hh,v 1.1 2003/08/26 18:41:50 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4Be9GEMChannel_h
#define G4Be9GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4Be9GEMCoulombBarrier.hh"
#include "G4Be9GEMProbability.hh"

class G4Be9GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4Be9GEMChannel() : G4GEMChannel(9,4,"Be9",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4Be9GEMChannel() {};
  
private:
  const G4Be9GEMChannel & operator=(const G4Be9GEMChannel & right);  
    
  G4Be9GEMChannel(const G4Be9GEMChannel & right);
  
public:
  G4bool operator==(const G4Be9GEMChannel & right) const;
  G4bool operator!=(const G4Be9GEMChannel & right) const;
    
private:
  
  G4Be9GEMCoulombBarrier theCoulombBarrier;
	
  G4Be9GEMProbability theEvaporationProbability;
  
};
#endif
