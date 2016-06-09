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
#include "G4DiffractiveHHScatterer.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4ExcitedString.hh"
#include "G4LundStringFragmentation.hh"
#include "G4KineticTrack.hh"
#include "G4DiffractiveSplitableHadron.hh"

G4DiffractiveHHScatterer::G4DiffractiveHHScatterer()
:
  theExcitation(new G4DiffractiveExcitation()),
  theStringFragmentation(new G4LundStringFragmentation())
{}

G4KineticTrackVector * G4DiffractiveHHScatterer::
Scatter(const G4KineticTrack & aTrack, const G4KineticTrack & bTrack)
{
  G4KineticTrackVector * result = new G4KineticTrackVector();

  G4DiffractiveSplitableHadron aHadron(& aTrack);
  G4DiffractiveSplitableHadron bHadron(& bTrack);
  if ( ! theExcitation->ExciteParticipants(& aHadron, & bHadron)) 
  {
	return NULL;
  }
  
  G4ExcitedString * string;
  G4KineticTrackVector * fragments=NULL;

  string = theExcitation->String(& aHadron, true);  // Projectile
  fragments=theStringFragmentation->FragmentString(*string);
  for (unsigned int aFragment=0; aFragment < fragments->size(); aFragment++)
  {
  	result->push_back(fragments->operator[](aFragment));
  }
  
  string = theExcitation->String(& bHadron, false); // Target
  fragments=theStringFragmentation->FragmentString(*string);
  for (unsigned int bFragment=0; bFragment < fragments->size(); bFragment++)
  {
  	result->push_back(fragments->operator[](bFragment));
  }
  
  
  return result;
}
