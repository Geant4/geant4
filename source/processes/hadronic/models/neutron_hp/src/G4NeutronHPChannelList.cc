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
    
  G4ParticleChange * G4NeutronHPChannelList::ApplyYourself(const G4Element * anElement, const G4Track & aTrack)
  {
    G4int i, ii;
    // decide on the isotope
    G4int numberOfIsos;
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
	  running[i] +=theChannels[ii]->GetWeightedXsec(aTrack.GetKineticEnergy(), i);
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
        running[i] += theChannels[i]->GetFSCrossSection(aTrack.GetKineticEnergy(), isotope);
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
//    G4cout << theDir << endl;
    theElement = anElement;
//    G4cout << theElement << endl;
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
