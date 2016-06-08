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
// $Id: G4C14GEMChannel.hh,v 1.1 2002/06/06 17:40:36 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4C14GEMChannel_h
#define G4C14GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4C14GEMCoulombBarrier.hh"
#include "G4C14GEMProbability.hh"

class G4C14GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4C14GEMChannel() : G4GEMChannel(14,6,"C14",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4C14GEMChannel() {};
  
private:
  const G4C14GEMChannel & operator=(const G4C14GEMChannel & right);  
    
  G4C14GEMChannel(const G4C14GEMChannel & right);
  
public:
  G4bool operator==(const G4C14GEMChannel & right) const;
  G4bool operator!=(const G4C14GEMChannel & right) const;
    
private:
  
  G4C14GEMCoulombBarrier theCoulombBarrier;
	
  G4C14GEMProbability theEvaporationProbability;
  
};
#endif
