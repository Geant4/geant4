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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPChannelList.hh"
#include "G4Element.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4NeutronHPFinalState.hh"

  G4int G4NeutronHPChannelList::trycounter = 0;

  G4NeutronHPChannelList::G4NeutronHPChannelList(G4int n)
  { 
    nChannels = n;
    theChannels  = new G4NeutronHPChannel * [n];
    allChannelsCreated = false;
    theInitCount = 0;
  }
  
  G4NeutronHPChannelList::G4NeutronHPChannelList()
  {
    nChannels = 0;
    theChannels = NULL;
    allChannelsCreated = false;
    theInitCount = 0;
  }
  
  G4NeutronHPChannelList::~G4NeutronHPChannelList()
  {
    if(theChannels!=NULL)
    {
      for(G4int i=0;i<nChannels; i++)
      {
        delete theChannels[i];
      }
      delete [] theChannels;
    }
  }
    
  #include "G4NeutronHPThermalBoost.hh"
  G4ParticleChange * G4NeutronHPChannelList::ApplyYourself(const G4Element * anElement, const G4Track & aTrack)
  {
    G4NeutronHPThermalBoost aThermalE;
    G4int i, ii;
    // decide on the isotope
    G4int numberOfIsos(0);
    for(ii=0; ii<nChannels; ii++)
    {
      numberOfIsos = theChannels[ii]->GetNiso();
      if(numberOfIsos!=0) break;
    }
    G4double * running= new G4double [numberOfIsos];
    running[0] = 0;
    for(i=0;i<numberOfIsos; i++)
    {
      if(i!=0) running[i] = running[i-1];
      for(ii=0; ii<nChannels; ii++)
      {
	if(theChannels[ii]->HasAnyData(i))
	{
          running[i] +=theChannels[ii]->GetWeightedXsec(aThermalE.GetThermalEnergy(aTrack.GetDynamicParticle(),
		                                                  theChannels[ii]->GetN(i),
								  theChannels[ii]->GetZ(i),
						  		  aTrack.GetMaterial()->GetTemperature()),
					                i);
	}
      }
    }
    G4int isotope=nChannels-1;
    G4double random=G4UniformRand();
    for(i=0;i<numberOfIsos; i++)
    {
      isotope = i;
      if(random<running[i]/running[numberOfIsos-1]) break;
    }
    delete [] running;
    
    // decide on the channel
    running = new G4double[nChannels];
    running[0]=0;
    for(i=0; i<nChannels; i++)
    {
      if(i!=0) running[i] = running[i-1];
      if(theChannels[i]->HasAnyData(isotope))
      {
        running[i] += theChannels[i]->GetFSCrossSection(aThermalE.GetThermalEnergy(aTrack.GetDynamicParticle(),
		                                                  theChannels[i]->GetN(isotope),
								  theChannels[i]->GetZ(isotope),
						  		  aTrack.GetMaterial()->GetTemperature()),
					                isotope);
      }
    }
    G4int lChan=0;
    random=G4UniformRand();
    for(i=0; i<nChannels; i++)
    {
      lChan = i;
      if(random<running[i]/running[nChannels-1]) break;
    }
    delete [] running;
    return theChannels[lChan]->ApplyYourself(aTrack, isotope);
  }
      
  void G4NeutronHPChannelList::Init(G4Element * anElement, const G4String & dirName)
  {
    theDir = dirName;
//    G4cout << theDir << G4endl;
    theElement = anElement;
//    G4cout << theElement << G4endl;
    ;
  }
  
  void G4NeutronHPChannelList::Register(G4NeutronHPFinalState * theFS, 
                                        const G4String & aName)
  {
    G4bool result;
    if(!allChannelsCreated)
    {
      if(nChannels!=0)
      {
	G4NeutronHPChannel ** theBuffer = new G4NeutronHPChannel * [nChannels+1];
	G4int i;
	for(i=0; i<nChannels; i++)
	{
	  theBuffer[i] = theChannels[i];
	}
	delete [] theChannels;
	theChannels = theBuffer;
      }
      else
      {
	theChannels = new G4NeutronHPChannel * [nChannels+1];
      }
      G4String name;
      name = aName+"/";
      theChannels[nChannels] = new G4NeutronHPChannel;
      theChannels[nChannels]->Init(theElement, theDir, name);
      nChannels++;
    }
    result = theChannels[theInitCount]->Register(theFS);
    theInitCount++; 
  }
