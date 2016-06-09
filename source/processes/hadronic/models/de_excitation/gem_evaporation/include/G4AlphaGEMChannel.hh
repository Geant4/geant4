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
// $Id: G4AlphaGEMChannel.hh,v 1.1 2003/08/26 18:41:35 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4AlphaGEMChannel_h
#define G4AlphaGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4AlphaCoulombBarrier.hh"
#include "G4AlphaGEMProbability.hh"

class G4AlphaGEMChannel : public G4GEMChannel
{
public:
    // only available constructor
    G4AlphaGEMChannel() : G4GEMChannel(4,2,"Alpha",
                                     &theEvaporationProbability,
                                     &theCoulombBarrier)
        {
            theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
        }

    // destructor
    ~G4AlphaGEMChannel() {};

private:
    const G4AlphaGEMChannel & operator=(const G4AlphaGEMChannel & right);  
    
    G4AlphaGEMChannel(const G4AlphaGEMChannel & right);
    
public:
    G4bool operator==(const G4AlphaGEMChannel & right) const;
    G4bool operator!=(const G4AlphaGEMChannel & right) const;
    
private:
    
    G4AlphaCoulombBarrier theCoulombBarrier;
	
    G4AlphaGEMProbability theEvaporationProbability;
    
};
#endif
