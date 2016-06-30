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
// $Id: G4Be7GEMChannel.hh 97097 2016-05-25 07:50:31Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept. 2001)
//
#ifndef G4Be7GEMChannel_h
#define G4Be7GEMChannel_h 1

#include "G4GEMChannel.hh"
#include "G4Be7GEMProbability.hh"

class G4Be7GEMChannel : public G4GEMChannel
{
public:
  // only available constructor
  explicit G4Be7GEMChannel() 
    : G4GEMChannel(7,4,"Be7",&theEvaporationProbability) {}
  
  // destructor
  ~G4Be7GEMChannel() {};
  
private:
  const G4Be7GEMChannel & operator=(const G4Be7GEMChannel & right) = delete;  
  G4Be7GEMChannel(const G4Be7GEMChannel & right) = delete;
  G4bool operator==(const G4Be7GEMChannel & right) const = delete;
  G4bool operator!=(const G4Be7GEMChannel & right) const = delete;
    
  G4Be7GEMProbability theEvaporationProbability;
  
};
#endif
