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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//
#ifndef G4C13GEMChannel_h
#define G4C13GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4C13GEMProbability.hh"

class G4C13GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  explicit G4C13GEMChannel() 
    : G4GEMChannel(13,6,"C13",&theEvaporationProbability) {}
  
  // destructor
  ~G4C13GEMChannel() {};
  
private:
  const G4C13GEMChannel & operator=(const G4C13GEMChannel & right) = delete;  
  G4C13GEMChannel(const G4C13GEMChannel & right) = delete;
  G4bool operator==(const G4C13GEMChannel & right) const = delete;
  G4bool operator!=(const G4C13GEMChannel & right) const = delete;
    
  G4C13GEMProbability theEvaporationProbability;
  
};
#endif
