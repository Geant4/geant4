// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPEnergyDistribution.hh,v 1.1 1999-01-07 16:12:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPEnergyDistribution_h
#define G4NeutronHPEnergyDistribution_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4NeutronHPArbitaryTab.hh"
#include "G4NeutronHPEvapSpectrum.hh"
#include "G4NeutronHPSimpleEvapSpectrum.hh"
#include "G4NeutronHPFissionSpectrum.hh"
#include "G4NeutronHPWattSpectrum.hh"
#include "G4NeutronHPMadlandNixSpectrum.hh"
#include "G4VNeutronHPEDis.hh"
#include "Randomize.hh"

// we will need a List of these .... one per term.

class G4NeutronHPEnergyDistribution
{
  public:
  G4NeutronHPEnergyDistribution()
  {
    theEnergyDistribution = NULL;
    theNumberOfPartials = 0;
    theRepresentationType = 0;
  }
  ~G4NeutronHPEnergyDistribution()
  {
    if(theEnergyDistribution != NULL)
    {
      for(G4int i=0; i<theNumberOfPartials; i++) 
      {
        delete theEnergyDistribution[i];
      }
      delete [] theEnergyDistribution;
    }
  }
  
  inline void Init(ifstream & theData)
  {
    G4double dummy;
    theData >> dummy >> theNumberOfPartials;
    theEnergyDistribution = new G4VNeutronHPEDis * [theNumberOfPartials];
    for(G4int i=0; i<theNumberOfPartials; i++) 
    {
      theData >> theRepresentationType;
      switch(theRepresentationType)
      {
	case 1:
          theEnergyDistribution[i] = new G4NeutronHPArbitaryTab;
          break;
	case 5:        
          theEnergyDistribution[i] = new G4NeutronHPEvapSpectrum;
          break;
	case 7:
          theEnergyDistribution[i] = new G4NeutronHPFissionSpectrum;
          break;
	case 9:
          theEnergyDistribution[i] = new G4NeutronHPSimpleEvapSpectrum;
          break;
	case 11:
          theEnergyDistribution[i] = new G4NeutronHPWattSpectrum;
          break;
	case 12:
          theEnergyDistribution[i] = new G4NeutronHPMadlandNixSpectrum;
          break;
      }
      theEnergyDistribution[i]->Init(theData);
    }
  }
  
  inline G4double Sample(G4double anEnergy, G4int & it) 
  {
    G4double result = 0;
    it = 0;
    if (theNumberOfPartials != 0)
    {
      G4double sum=0;
      G4double * running = new G4double[theNumberOfPartials];
      running[0] = 0;
      G4int i;
      for (i=0; i<theNumberOfPartials; i++)
      {	
	if (i!=0) running[i]=running[i-1];
	running[i]+=theEnergyDistribution[i]->GetFractionalProbability(anEnergy);
      }
      sum = running[theNumberOfPartials-1];
      G4double random = G4UniformRand();
      for(i=0; i<theNumberOfPartials; i++)
      {
	it = i;
	if(running[i]/sum>random) break;
      }
      delete [] running;
      if(it==theNumberOfPartials) it--;
      result = theEnergyDistribution[it]->Sample(anEnergy);
    }
    return result;
  }
  
  private:
  
  G4int theNumberOfPartials;
  G4int theRepresentationType;
  G4VNeutronHPEDis ** theEnergyDistribution;
};

#endif
