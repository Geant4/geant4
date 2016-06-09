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
// $Id: G4B10GEMChannel.hh,v 1.2 2005/06/04 13:25:24 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4B10GEMChannel_h
#define G4B10GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4B10GEMCoulombBarrier.hh"
#include "G4B10GEMProbability.hh"

class G4B10GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4B10GEMChannel() : G4GEMChannel(10,5,"B10",
				   &theEvaporationProbability,
				   &theCoulombBarrier)
  {
    theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
  }
  
  // destructor
  ~G4B10GEMChannel() {};
  
private:
  const G4B10GEMChannel & operator=(const G4B10GEMChannel & right);  
    
  G4B10GEMChannel(const G4B10GEMChannel & right);
  
public:
  G4bool operator==(const G4B10GEMChannel & right) const;
  G4bool operator!=(const G4B10GEMChannel & right) const;
    
private:
  
  G4B10GEMCoulombBarrier theCoulombBarrier;
	
  G4B10GEMProbability theEvaporationProbability;
  
};
#endif
