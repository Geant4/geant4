//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4DeuteronGEMChannel.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//
// J. M. Quesada (July 2009) coulomb barrier striclty according to Furihata's paper

#ifndef G4DeuteronGEMChannel_h
#define G4DeuteronGEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4DeuteronGEMCoulombBarrier.hh"
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
     // JMQ 190709
//     G4DeuteronCoulombBarrier theCoulombBarrier;    
    G4DeuteronGEMCoulombBarrier theCoulombBarrier;
	
    G4DeuteronGEMProbability theEvaporationProbability;
    
};
#endif
