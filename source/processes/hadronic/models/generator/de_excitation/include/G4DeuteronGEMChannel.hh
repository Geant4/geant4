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
// $Id: G4DeuteronGEMChannel.hh,v 1.1 2002/06/06 17:41:07 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4DeuteronGEMChannel_h
#define G4DeuteronGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4DeuteronCoulombBarrier.hh"
#include "G4DeuteronGEMProbability.hh"

class G4DeuteronGEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4DeuteronGEMChannel() : G4GEMChannel(2,1,"deuteron",
                                        &theEvaporationProbability,
                                        &theCoulombBarrier)
        {
            theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
            
        }
  // destructor
  ~G4DeuteronGEMChannel() {};

private:
    const G4DeuteronGEMChannel & operator=(const G4DeuteronGEMChannel & right);  
    
    G4DeuteronGEMChannel(const G4DeuteronGEMChannel & right);
    
public:
    G4bool operator==(const G4DeuteronGEMChannel & right) const;
    G4bool operator!=(const G4DeuteronGEMChannel & right) const;
    
private:
    
    G4DeuteronCoulombBarrier theCoulombBarrier;
	
    G4DeuteronGEMProbability theEvaporationProbability;
    
};
#endif
