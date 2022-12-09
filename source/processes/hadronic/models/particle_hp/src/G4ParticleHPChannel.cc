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
// 071031 bug fix T. Koi on behalf of A. Howard 
// 081203 bug fix in Register method by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
// June-2019 - E. Mendoza --> Modification to allow using an incomplete
//   data library if the G4NEUTRONHP_SKIP_MISSING_ISOTOPES environmental
//   flag is defined. The missing XS are set to 0.


#include <stdlib.h>

#include "G4ParticleHPChannel.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4HadTmpUtil.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4ParticleHPReactionWhiteBoard.hh"

G4double G4ParticleHPChannel::GetXsec(G4double energy)
{
  return std::max(0., theChannelData->GetXsec(energy));
}
  
G4double G4ParticleHPChannel::GetWeightedXsec(G4double energy, G4int isoNumber)
{
  return theIsotopeWiseData[isoNumber].GetXsec(energy);
}
  
G4double G4ParticleHPChannel::GetFSCrossSection(G4double energy, G4int isoNumber)
{
  return theFinalStates[isoNumber]->GetXsec(energy);
}
  
  void G4ParticleHPChannel::
  Init(G4Element * anElement, const G4String dirName, const G4String aFSType) 
  {
    theFSType = aFSType;
    Init(anElement, dirName);
  }
   
  void G4ParticleHPChannel::Init(G4Element * anElement, const G4String dirName)  
  {
    theDir = dirName;
    theElement = anElement;
  }
  
  G4bool G4ParticleHPChannel::Register(G4ParticleHPFinalState *theFS)
  {
	registerCount++;
    G4int Z = G4lrint(theElement->GetZ());

    Z = Z-registerCount;
    if ( registerCount > 5 ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this material"); // for Elastic, Capture, Fission case 
    if ( Z < 1 ) return false; 
/*
    if(registerCount<5)
    {
      Z = Z-registerCount;
    }
*/
    //if(Z=theElement->GetZ()-5) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this material");
    // Bug fix by TK on behalf of AH
    //if ( Z <=theElement->GetZ()-5 ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this material");
    G4int count = 0;
    if(registerCount==0) count = (G4int)theElement->GetNumberOfIsotopes();
    if(count == 0||registerCount!=0) count +=
         theStableOnes.GetNumberOfIsotopes(Z);
    niso = count;
    delete [] theIsotopeWiseData;
    theIsotopeWiseData = new G4ParticleHPIsoData [niso];
    delete [] active;
    active = new G4bool[niso];

    delete [] theFinalStates;
    theFinalStates = new G4ParticleHPFinalState * [niso];
    delete theChannelData;
    theChannelData = new G4ParticleHPVector; 
    for(G4int i=0; i<niso; ++i)
    {
      theFinalStates[i] = theFS->New();
      theFinalStates[i]->SetProjectile(theProjectile);
    }
    count = 0;
    G4int nIsos = niso;
    if(theElement->GetNumberOfIsotopes()!=0&&registerCount==0)
    {
      for (G4int i1=0; i1<nIsos; ++i1)
      {
        // G4cout <<" Init: normal case"<<G4endl;
        G4int A = theElement->GetIsotope(i1)->GetN();
        G4int M = theElement->GetIsotope(i1)->Getm();
        G4double frac = theElement->GetRelativeAbundanceVector()[i1]/perCent;
        //theFinalStates[i1]->SetA_Z(A, Z);
	//UpdateData(A, Z, count++, frac);
        theFinalStates[i1]->SetA_Z(A, Z, M);
	UpdateData(A, Z, M, count++, frac, theProjectile);
      }
    } else {
      //G4cout <<" Init: mean case: "
      //       <<theStableOnes.GetNumberOfIsotopes(Z)<<" "
	//     <<Z<<" "<<theElement
	//     << G4endl;
      G4int first = theStableOnes.GetFirstIsotope(Z);
      for(G4int i1=0; i1<theStableOnes.GetNumberOfIsotopes(Z); ++i1)
      {
        G4int A = theStableOnes.GetIsotopeNucleonCount(first+i1);
        G4double frac = theStableOnes.GetAbundance(first+i1);
        theFinalStates[i1]->SetA_Z(A, Z);
        UpdateData(A, Z, count++, frac, theProjectile);
      }
    }
    G4bool result = HasDataInAnyFinalState();

    //To avoid issuing hash by worker threads
    if ( result ) theChannelData->Hash();

    return result;
  }
  
  void G4ParticleHPChannel::UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance, G4ParticleDefinition* projectile)
  {
    // Initialze the G4FissionFragment generator for this isomer if needed
    if(wendtFissionGenerator)
    {
      wendtFissionGenerator->InitializeANucleus(A, Z, M, theDir);
    }

    theFinalStates[index]->Init(A, Z, M, theDir, theFSType, projectile);
    if(!theFinalStates[index]->HasAnyData()) return; // nothing there for exactly this isotope.

    // the above has put the X-sec into the FS
    theBuffer = 0;
    if(theFinalStates[index]->HasXsec())
    {
      theBuffer = theFinalStates[index]->GetXsec();
      theBuffer->Times(abundance/100.);
      theIsotopeWiseData[index].FillChannelData(theBuffer);
    }
    else // get data from CrossSection directory
    {
      G4String tString = "/CrossSection";
      //active[index] = theIsotopeWiseData[index].Init(A, Z, abundance, theDir, tString);
      active[index] = theIsotopeWiseData[index].Init(A, Z, M, abundance, theDir, tString);
      if(active[index]) theBuffer = theIsotopeWiseData[index].MakeChannelData();
    }
    if(theBuffer != 0) Harmonise(theChannelData, theBuffer);
  }
  
  void G4ParticleHPChannel::Harmonise(G4ParticleHPVector *& theStore, G4ParticleHPVector * theNew)
  {
    G4int s_tmp = 0, n=0, m_tmp=0;
    G4ParticleHPVector * theMerge = new G4ParticleHPVector;
    G4ParticleHPVector * anActive = theStore;
    G4ParticleHPVector * aPassive = theNew;
    G4ParticleHPVector * tmp;
    G4int a = s_tmp, p = n, t;
    while (a<anActive->GetVectorLength()&&p<aPassive->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      if(anActive->GetEnergy(a) <= aPassive->GetEnergy(p))
      {
        G4double xa  = anActive->GetEnergy(a);
        theMerge->SetData(m_tmp, xa, anActive->GetXsec(a)+std::max(0., aPassive->GetXsec(xa)) );
        m_tmp++;
        a++;
        G4double xp = aPassive->GetEnergy(p);
        if( std::abs(std::abs(xp-xa)/xa)<0.001 )
        {
          ++p;
        }
      } else {
        tmp = anActive; t=a;
        anActive = aPassive; a=p;
        aPassive = tmp; p=t;
      }
    }
    while (a!=anActive->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      theMerge->SetData(m_tmp++, anActive->GetEnergy(a), anActive->GetXsec(a));
      ++a;
    }
    while (p!=aPassive->GetVectorLength()) // Loop checking, 11.05.2015, T. Koi
    {
      if(std::abs(theMerge->GetEnergy(std::max(0,m_tmp-1))-aPassive->GetEnergy(p))/aPassive->GetEnergy(p)>0.001)
        theMerge->SetData(m_tmp++, aPassive->GetEnergy(p), aPassive->GetXsec(p));
      ++p;
    }
    delete theStore;
    theStore = theMerge;
  }

  G4HadFinalState * G4ParticleHPChannel::
  ApplyYourself(const G4HadProjectile & theTrack, G4int anIsotope)
  {
//    G4cout << "G4ParticleHPChannel::ApplyYourself+"<<niso<<G4endl;
    if ( anIsotope != -1 && anIsotope != -2 ) 
    {
       //Inelastic Case
       //G4cout << "G4ParticleHPChannel Inelastic Case" 
       //<< " Z= " << this->GetZ(it) << " A = " << this->GetN(it) << G4endl;
       G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargA( (G4int)this->GetN(anIsotope) ); 
       G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargZ( (G4int)this->GetZ(anIsotope) ); 
       return theFinalStates[anIsotope]->ApplyYourself(theTrack);
    }
    G4double sum=0;
    G4int it=0;
    G4double * xsec = new G4double[niso];
    G4ParticleHPThermalBoost aThermalE;
    for (G4int i=0; i<niso; i++)
    {
      if(theFinalStates[i]->HasAnyData())
      {
        xsec[i] = theIsotopeWiseData[i].GetXsec(aThermalE.GetThermalEnergy(theTrack,
		                                                           theFinalStates[i]->GetN(),
									   theFinalStates[i]->GetZ(),
						  		           theTrack.GetMaterial()->GetTemperature()));
        sum += xsec[i];
      }
      else
      {
        xsec[i]=0;
      }
    } 
    if(sum == 0) 
    {
//      G4cout << "G4ParticleHPChannel::ApplyYourself theFinalState->Initialize+"<<G4endl;
//      G4cout << "G4ParticleHPChannel::ApplyYourself theFinalState->Initialize-"<<G4endl;
      it = static_cast<G4int>(niso*G4UniformRand());
    }
    else
    {
//      G4cout << "Are we still here? "<<sum<<G4endl;
//      G4cout << "TESTHP 23 NISO="<<niso<<G4endl;
      G4double random = G4UniformRand();
      G4double running=0;
//      G4cout << "G4ParticleHPChannel::ApplyYourself Done the sum"<<niso<<G4endl;
//      G4cout << "TESTHP 24 NISO="<<niso<<G4endl;
      for (G4int ix=0; ix<niso; ix++)
      {
        running += xsec[ix];
        //if(random<=running/sum) 
        if( sum == 0 || random <= running/sum ) 
        { 
          it = ix;
	  break;
        }
      }
      if(it==niso) it--;
    }
    delete [] xsec;
    G4HadFinalState * theFinalState=0;
    const G4int A = (G4int)this->GetN(it);
    const G4int Z = (G4int)this->GetZ(it);
    const G4int M = (G4int)this->GetM(it);

                                       //-2:Marker for Fission
    if(wendtFissionGenerator&&anIsotope==-2)
    {
    	theFinalState = wendtFissionGenerator->ApplyYourself(theTrack, Z, A);
    }

    // Use the standard procedure if the G4FissionFragmentGenerator model fails
    if (!theFinalState)
    {

       G4int icounter=0;
       G4int icounter_max=1024;
       while(theFinalState==0) // Loop checking, 11.05.2015, T. Koi
       {
          icounter++;
          if ( icounter > icounter_max ) {
	     G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
             break;
          }
//	      G4cout << "TESTHP 24 it="<<it<<G4endl;
          theFinalState = theFinalStates[it]->ApplyYourself(theTrack);
       }
    }

    //G4cout <<"THE IMPORTANT RETURN"<<G4endl;
    //G4cout << "TK G4ParticleHPChannel Elastic, Capture and Fission Cases "
    //<< " Z= " << this->GetZ(it) << " A = " << this->GetN(it) << G4endl;
    G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargA( A );
    G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargZ( Z );
    G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->SetTargM( M );

    return theFinalState;
  }


void G4ParticleHPChannel::DumpInfo(){

  G4cout<<" Element: "<<theElement->GetName()<<G4endl;
  G4cout<<" Directory name: "<<theDir<<G4endl;
  G4cout<<" FS name: "<<theFSType<<G4endl;
  G4cout<<" Number of Isotopes: "<<niso<<G4endl;
  G4cout<<" Have cross sections: "<<G4endl;
  for(int i=0;i<niso;i++){
    G4cout<<theFinalStates[i]->HasXsec()<<"  ";
  }
  G4cout<<G4endl;
  if(theChannelData){
    G4cout<<" Cross Section (total for this channel):"<<G4endl;
    int np=theChannelData->GetVectorLength();
    G4cout<<np<<G4endl;
    for(int i=0;i<np;i++){
      G4cout<<theChannelData->GetEnergy(i)/eV<<"  "<<theChannelData->GetXsec(i)<<G4endl;
    }
  }

}
 


