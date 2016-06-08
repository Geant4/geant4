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
// $Id: G4Evaporation.hh,v 1.5.2.1 2001/06/28 19:13:01 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"

#include "G4ios.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"

//#define debug

class G4Evaporation : public G4VEvaporation
{
public:
  G4Evaporation();
  G4Evaporation(G4std::vector<G4VEvaporationChannel*> * aChannelsVector) :
    theChannels(aChannelsVector),
    myOwnChannelsVector(false)
  {};
	 
  ~G4Evaporation();

private:
  G4Evaporation(const G4Evaporation &right);

  const G4Evaporation & operator=(const G4Evaporation &right);
  G4bool operator==(const G4Evaporation &right) const;
  G4bool operator!=(const G4Evaporation &right) const;

public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
		
private:

#ifdef debug
  void CheckConservation(const G4Fragment & theInitialState,
			 G4FragmentVector * Result) const;
#endif

  G4bool myOwnChannelsVector;

  G4std::vector<G4VEvaporationChannel*> * theChannels;
};

#endif





