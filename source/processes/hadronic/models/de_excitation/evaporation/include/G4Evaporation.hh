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
// $Id: G4Evaporation.hh,v 1.12 2010-05-11 11:34:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
//
// 14/02/2007 Alex Howard - added protection for negative probabilities in the sum 
// 27/07/2009 V.Ivanchenko - added Combined decay channels (default + GEM) 
// 11/05/2010 V.Ivanchenko - rewrited technical part do not "new" and "delete" 
//                           of small objects

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4UnstableFragmentBreakUp.hh"

class G4VEvaporationFactory;
class G4NistManager;

//#define debug

class G4Evaporation : public G4VEvaporation
{
public:
  G4Evaporation();
  G4Evaporation(std::vector<G4VEvaporationChannel*> * aChannelsVector);
	 
  virtual ~G4Evaporation();

  virtual void Initialise();

private:

  void InitialiseEvaporation();

  G4Evaporation(const G4Evaporation &right);

  const G4Evaporation & operator=(const G4Evaporation &right);
  G4bool operator==(const G4Evaporation &right) const;
  G4bool operator!=(const G4Evaporation &right) const;

public:

  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);

  void SetDefaultChannel();
  void SetGEMChannel();
  void SetCombinedChannel();

#ifdef debug
  void CheckConservation(const G4Fragment & theInitialState,
			 G4FragmentVector * Result) const;
#endif

  std::vector<G4VEvaporationChannel*> * theChannels;
  std::vector<G4double>   probabilities;
  G4VEvaporationFactory * theChannelFactory;
  G4int nChannels;
  G4double minExcitation;
  G4NistManager* nist;
  G4UnstableFragmentBreakUp unstableBreakUp;
    
};

#endif





