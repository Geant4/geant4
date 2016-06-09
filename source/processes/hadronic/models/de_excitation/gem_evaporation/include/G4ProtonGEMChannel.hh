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
// $Id: G4ProtonGEMChannel.hh,v 1.1 2003/08/26 18:43:16 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//


#ifndef G4ProtonGEMChannel_h
#define G4ProtonGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4ProtonCoulombBarrier.hh"
#include "G4ProtonGEMProbability.hh"

class G4ProtonGEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  G4ProtonGEMChannel() : G4GEMChannel(1,1,"proton",
                                      &theEvaporationProbability,
                                      &theCoulombBarrier)
        {
            theEvaporationProbability.SetCoulomBarrier(&theCoulombBarrier);
        }
    
  // destructor
  ~G4ProtonGEMChannel() {};

private:
  const G4ProtonGEMChannel & operator=(const G4ProtonGEMChannel & right);  

  G4ProtonGEMChannel(const G4ProtonGEMChannel & right);

public:
  G4bool operator==(const G4ProtonGEMChannel & right) const;
  G4bool operator!=(const G4ProtonGEMChannel & right) const;

private:

  G4ProtonCoulombBarrier theCoulombBarrier;
	
  G4ProtonGEMProbability theEvaporationProbability;

};
#endif
