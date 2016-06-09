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
// $Id: G4TritonGEMChannel.hh,v 1.1 2002/06/06 17:53:47 larazb Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4TritonGEMChannel_h
#define G4TritonGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4TritonCoulombBarrier.hh"
#include "G4TritonGEMProbability.hh"

class G4TritonGEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4TritonGEMChannel() : G4GEMChannel(3,1,"triton",
                                      &theEvaporationProbability,
                                      &theCoulombBarrier)
        {
            theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
        }
    
  // destructor
  ~G4TritonGEMChannel() {};

private:
    const G4TritonGEMChannel & operator=(const G4TritonGEMChannel & right);  
    
    G4TritonGEMChannel(const G4TritonGEMChannel & right);
    
public:
    G4bool operator==(const G4TritonGEMChannel & right) const;
    G4bool operator!=(const G4TritonGEMChannel & right) const;
    
private:
    
    G4TritonCoulombBarrier theCoulombBarrier;
	
    G4TritonGEMProbability theEvaporationProbability;
    
};
#endif
