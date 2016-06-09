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
// $Id: G4AlphaEvaporationChannel.hh,v 1.3 2006/06/29 20:09:43 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov. 1999)
//


#ifndef G4AlphaEvaporationChannel_h
#define G4AlphaEvaporationChannel_h 1

#include "G4EvaporationChannel.hh"
#include "G4AlphaCoulombBarrier.hh"
#include "G4AlphaEvaporationProbability.hh"

class G4AlphaEvaporationChannel : public G4EvaporationChannel
{
public:
  // only available constructor
  G4AlphaEvaporationChannel() : G4EvaporationChannel(4,2,"alpha",
						     &theEvaporationProbability,&theCoulombBarrier) {};

  // destructor
  ~G4AlphaEvaporationChannel() {};

private:
  const G4AlphaEvaporationChannel & operator=(const G4AlphaEvaporationChannel & right);  

  G4AlphaEvaporationChannel(const G4AlphaEvaporationChannel & right);


public:
  G4bool operator==(const G4AlphaEvaporationChannel & right) const;
  G4bool operator!=(const G4AlphaEvaporationChannel & right) const;

private:

  G4AlphaCoulombBarrier theCoulombBarrier;
	
  G4AlphaEvaporationProbability theEvaporationProbability;

};
#endif
