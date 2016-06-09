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
// $Id: G4Evaporation.hh,v 1.5 2008/09/19 13:32:54 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
//
// Alex Howard - added protection for negative probabilities in the sum, 14/2/07

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"

class G4VEvaporationFactory;

//#define debug

class G4Evaporation : public G4VEvaporation
{
public:
  G4Evaporation();
  G4Evaporation(std::vector<G4VEvaporationChannel*> * aChannelsVector) :
    theChannels(aChannelsVector), theChannelFactory(0)
  {};
	 
  ~G4Evaporation();

private:
  G4Evaporation(const G4Evaporation &right);

  const G4Evaporation & operator=(const G4Evaporation &right);
  G4bool operator==(const G4Evaporation &right) const;
  G4bool operator!=(const G4Evaporation &right) const;

public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);

  void SetDefaultChannel();
  void SetGEMChannel();

#ifdef debug
  void CheckConservation(const G4Fragment & theInitialState,
			 G4FragmentVector * Result) const;
#endif


  std::vector<G4VEvaporationChannel*> * theChannels;
  G4VEvaporationFactory * theChannelFactory;
  
  
  class SumProbabilities : public std::binary_function<G4double,G4double,G4double>
  {
  public:
    SumProbabilities() : total(0.0) {}
    G4double operator() (G4double& /* probSoFar */, G4VEvaporationChannel*& frag)
    { 
      G4double temp_prob = frag->GetEmissionProbability();
      if(temp_prob >= 0.0) total += temp_prob;
      //      total += frag->GetEmissionProbability();
      return total;
    }
    
    G4double GetTotal() { return total; }
  public:
    G4double total;
    
  };

};

#endif





