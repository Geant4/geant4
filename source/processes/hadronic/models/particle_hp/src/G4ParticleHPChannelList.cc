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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPChannelList.hh"
#include "G4Element.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4ParticleHPFinalState.hh"

G4ThreadLocal G4int G4ParticleHPChannelList::trycounter = 0;

 G4ParticleHPChannelList::G4ParticleHPChannelList(G4int n, G4ParticleDefinition* projectile)
  :theProjectile(projectile)
  { 
    nChannels = n;
    theChannels  = new G4ParticleHPChannel * [n];
    allChannelsCreated = false;
    theInitCount = 0;
    theElement = NULL;
  }
  
#include "G4Neutron.hh"
  G4ParticleHPChannelList::G4ParticleHPChannelList()
  {
    nChannels = 0;
    theChannels = 0;
    allChannelsCreated = false;
    theInitCount = 0;
    theElement = NULL;
    theProjectile = G4Neutron::Neutron();
  }
  
  G4ParticleHPChannelList::~G4ParticleHPChannelList()
  {
    if(theChannels!=0)
    {
      for(G4int i=0;i<nChannels; i++)
      {
        delete theChannels[i];
      }
      delete [] theChannels;
    }
  }
    
  #include "G4ParticleHPThermalBoost.hh"
  #include "G4ParticleHPManager.hh"
  G4HadFinalState * G4ParticleHPChannelList::ApplyYourself(const G4Element * , const G4HadProjectile & aTrack)
  {
    G4ParticleHPThermalBoost aThermalE;
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
          running[i] +=theChannels[ii]->GetWeightedXsec(aThermalE.GetThermalEnergy(aTrack,
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
      //if(random<running[i]/running[numberOfIsos-1]) break;
      if(running[numberOfIsos-1] != 0) if(random<running[i]/running[numberOfIsos-1]) break;
    }
    delete [] running;
    
     // decide on the channel
    running = new G4double[nChannels];
    running[0]=0;
    G4int targA=-1; // For production of unChanged
    G4int targZ=-1;
    for(i=0; i<nChannels; i++)
    {
      if(i!=0) running[i] = running[i-1];
      if(theChannels[i]->HasAnyData(isotope))
      {
        targA=(G4int)theChannels[i]->GetN(isotope); //Will be simply used the last valid value
        targZ=(G4int)theChannels[i]->GetZ(isotope);
        running[i] += theChannels[i]->GetFSCrossSection(aThermalE.GetThermalEnergy(aTrack,
		                                                  targA,
								  targZ,
						  		  aTrack.GetMaterial()->GetTemperature()),
					                isotope);
        targA=(G4int)theChannels[i]->GetN(isotope); //Will be simply used the last valid value
        targZ=(G4int)theChannels[i]->GetZ(isotope);
	//	G4cout << " G4ParticleHPChannelList " << i << " targA " << targA << " targZ " << targZ << " running " << running[i] << G4endl;
      }
    }

    //TK120607
    if ( running[nChannels-1] == 0 )
    {
       //This happened usually by the miss match between the cross section data and model
       if ( targA == -1 && targZ == -1 ) {
          throw G4HadronicException(__FILE__, __LINE__, "ParticleHP model encounter lethal discrepancy with cross section data");
       }

       //TK121106
       G4cout << "Warning from NeutronHP: could not find proper reaction channel. This may cause by inconsistency between cross section and model.  Unchanged final states are returned." << G4endl;
       unChanged.Clear();

       //For Ep Check create unchanged final state including rest target 
       G4ParticleDefinition* targ_pd = G4IonTable::GetIonTable()->GetIon ( targZ , targA , 0.0 );
       G4DynamicParticle* targ_dp = new G4DynamicParticle( targ_pd , G4ThreeVector(1,0,0), 0.0 );
       unChanged.SetEnergyChange(aTrack.GetKineticEnergy());
       unChanged.SetMomentumChange(aTrack.Get4Momentum().vect() );
       unChanged.AddSecondary(targ_dp);
       //TK121106
       G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargA( targA ); 
       G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargZ( targZ ); 
       delete [] running;
       return &unChanged;
    }
    //TK120607


    G4int lChan=0;
    random=G4UniformRand();
    for(i=0; i<nChannels; i++)
    {
      lChan = i;
      if(running[nChannels-1] != 0) if(random<running[i]/running[nChannels-1]) break;
    }
    delete [] running;
#ifdef G4PHPDEBUG
    if( getenv("G4ParticleHPDebug") ) G4cout << " G4ParticleHPChannelList SELECTED ISOTOPE " << isotope << " SELECTED CHANNEL " << lChan << G4endl;
#endif
    return theChannels[lChan]->ApplyYourself(aTrack, isotope);
  }
      
void G4ParticleHPChannelList::Init(G4Element * anElement, const G4String & dirName, G4ParticleDefinition* projectile )
  {
    theDir = dirName;
//    G4cout << theDir << G4endl;
    theElement = anElement;
//    G4cout << theElement << G4endl;
    theProjectile = projectile;
  }
  
  void G4ParticleHPChannelList::Register(G4ParticleHPFinalState * theFS, 
					 const G4String & aName )
  {
    if(!allChannelsCreated)
    {
      if(nChannels!=0)
      {
	G4ParticleHPChannel ** theBuffer = new G4ParticleHPChannel * [nChannels+1];
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
	theChannels = new G4ParticleHPChannel * [nChannels+1];
      }
      G4String name;
      name = aName+"/";
      theChannels[nChannels] = new G4ParticleHPChannel(theProjectile);
      theChannels[nChannels]->Init(theElement, theDir, name);
      //      theChannels[nChannels]->SetProjectile( theProjectile );
      nChannels++;
    }
    
    //110527TKDB  Unnessary codes, Detected by gcc4.6 compiler 
    //G4bool result;
    //result = theChannels[theInitCount]->Register(theFS);
    theChannels[theInitCount]->Register(theFS);
    theInitCount++; 
  }

void G4ParticleHPChannelList::DumpInfo(){

  G4cout<<"================================================================"<<G4endl;
  G4cout<<" Element: "<<theElement->GetName()<<G4endl;
  G4cout<<" Number of channels: "<<nChannels<<G4endl;
  G4cout<<" Projectile: "<<theProjectile->GetParticleName()<<G4endl;
  G4cout<<" Directory name: "<<theDir<<G4endl;
  for(int i=0;i<nChannels;i++){
    if(theChannels[i]->HasDataInAnyFinalState()){
      G4cout<<"----------------------------------------------------------------"<<G4endl;
      theChannels[i]->DumpInfo();
      G4cout<<"----------------------------------------------------------------"<<G4endl;
    }
  }
  G4cout<<"================================================================"<<G4endl;

}
