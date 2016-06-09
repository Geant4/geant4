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
// $Id: G4He3GEMChannel.hh,v 1.1 2003/08/26 18:42:13 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4He3GEMChannel_h
#define G4He3GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4He3CoulombBarrier.hh"
#include "G4He3GEMProbability.hh"

class G4He3GEMChannel : public G4GEMChannel
{
public:
    // only available constructor
    G4He3GEMChannel() : G4GEMChannel(3,2,"He3",
                                     &theEvaporationProbability,
                                     &theCoulombBarrier)
        {
            theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
        }

    // destructor
  ~G4He3GEMChannel() {};

private:
    const G4He3GEMChannel & operator=(const G4He3GEMChannel & right);  
    
    G4He3GEMChannel(const G4He3GEMChannel & right);
    
public:
    G4bool operator==(const G4He3GEMChannel & right) const;
    G4bool operator!=(const G4He3GEMChannel & right) const;
    
private:
    
    G4He3CoulombBarrier theCoulombBarrier;
	
    G4He3GEMProbability theEvaporationProbability;
    
};
#endif
