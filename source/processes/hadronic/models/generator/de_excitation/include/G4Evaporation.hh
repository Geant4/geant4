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
// $Id: G4Evaporation.hh,v 1.12 2002/12/12 19:17:05 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

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
  G4Evaporation(G4std::vector<G4VEvaporationChannel*> * aChannelsVector) :
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
  
private:

#ifdef debug
  void CheckConservation(const G4Fragment & theInitialState,
			 G4FragmentVector * Result) const;
#endif


  G4std::vector<G4VEvaporationChannel*> * theChannels;
  G4VEvaporationFactory * theChannelFactory;
  
  
  class SumProbabilities : public G4std::binary_function<G4double,G4double,G4double>
  {
  public:
    SumProbabilities() : total(0.0) {}
    G4double operator() (G4double& probSoFar, G4VEvaporationChannel*& frag)
    { 
      total += frag->GetEmissionProbability();
      return total;
    }
    
    G4double GetTotal() { return total; }
  public:
    G4double total;
    
  };

};

#endif





