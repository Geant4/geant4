// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPChannel.hh"
#include "G4NeutronHPFinalState.hh"
#include "globals.hh"

  G4double G4NeutronHPChannel::GetXsec(G4double energy)
  {
    return G4std::max(0., theChannelData->GetXsec(energy));
  }
  
  G4double G4NeutronHPChannel::GetWeightedXsec(G4double energy, G4int isoNumber)
  {
    return theIsotopeWiseData[isoNumber].GetXsec(energy);
  }
  
  G4double G4NeutronHPChannel::GetFSCrossSection(G4double energy, G4int isoNumber)
  {
    return theFinalStates[isoNumber]->GetXsec(energy);
  }
  
  void G4NeutronHPChannel::
  Init(G4Element * anElement, const G4String dirName, const G4String aFSType) 
  {
    theFSType = aFSType;
    Init(anElement, dirName);
  }
   
  void G4NeutronHPChannel::Init(G4Element * anElement, const G4String dirName)  
  {
    theDir = dirName;
    theElement = anElement;
  }
  
  G4bool G4NeutronHPChannel::Register(G4NeutronHPFinalState *theFS)
  {
    registerCount++;
    G4int Z = theElement->GetZ();
    if(registerCount<5)
    {
      Z = Z-registerCount;
    }
    if(Z==theElement->GetZ()-5) G4Exception("Channel: Do not know what to do with this material");
    G4int count = 0;
    if(registerCount==0) count = theElement->GetNumberOfIsotopes();
    if(count == 0||registerCount!=0) count +=
         theStableOnes.GetNumberOfIsotopes(Z);
    niso = count;
    if(theIsotopeWiseData!=NULL) delete [] theIsotopeWiseData;
    theIsotopeWiseData = new G4NeutronHPIsoData [niso];
    if(active!=NULL) delete [] active;
    active = new G4bool[niso];
    if(theFinalStates!=NULL) delete [] theFinalStates;
    theFinalStates = new G4NeutronHPFinalState * [niso];
    delete theChannelData;
    theChannelData = new G4NeutronHPVector; 
    for(G4int i=0; i<niso; i++)
    {
      theFinalStates[i] = theFS->New();
    }
    count = 0;
    G4int nIsos = niso;
    if(theElement->GetNumberOfIsotopes()!=0&&registerCount==0)
    {
      for (G4int i1=0; i1<nIsos; i1++)
      {
//        G4cout <<" Init: normal case"<<G4endl;
        G4int A = theElement->GetIsotope(i1)->GetN();
        G4double frac = theElement->GetRelativeAbundanceVector()[i1]/perCent;
        UpdateData(A, Z, count++, frac);
      }
    } else {
      G4int first = theStableOnes.GetFirstIsotope(Z);
      for(G4int i1=0; 
        i1<theStableOnes.GetNumberOfIsotopes(Z);
        i1++)
      {
        G4int A = theStableOnes.GetIsotopeNucleonCount(first+i1);
        G4double frac = theStableOnes.GetAbundance(first+i1);
        UpdateData(A, Z, count++, frac);
      }
    }
    G4bool result = HasDataInAnyFinalState();
    return result;
  }
  
  void G4NeutronHPChannel::UpdateData(G4int A, G4int Z, G4int index, G4double abundance)
  {
    theFinalStates[index]->Init(A, Z, theDir, theFSType);
    if(!theFinalStates[index]->HasAnyData()) return; // nothing there for exactly this isotope.

    // the above has put the X-sec into the FS
    theBuffer = NULL;
    if(theFinalStates[index]->HasXsec())
    {
      theBuffer = theFinalStates[index]->GetXsec();
      theBuffer->Times(abundance/100.);
      theIsotopeWiseData[index].FillChannelData(theBuffer);
    }
    else // get data from CrossSection directory
    {
      G4String tString = "/CrossSection/";
      active[index] = theIsotopeWiseData[index].Init(A, Z, abundance, theDir, tString);
      if(active[index]) theBuffer = theIsotopeWiseData[index].MakeChannelData();
    }
    if(theBuffer != NULL) Harmonise(theChannelData, theBuffer);
  }
  
  void G4NeutronHPChannel::Harmonise(G4NeutronHPVector *& theStore, G4NeutronHPVector * theNew)
  {
    G4int s = 0, n=0, i=0, m=0;
    G4NeutronHPVector * theMerge = new G4NeutronHPVector;
    G4bool flag;
    G4NeutronHPVector * anActive = theStore;
    G4NeutronHPVector * aPassive = theNew;
    G4NeutronHPVector * tmp;
    G4int a = s, p = n, t;
    while (a<anActive->GetVectorLength()&&p<aPassive->GetVectorLength())
    {
      if(anActive->GetEnergy(a) <= aPassive->GetEnergy(p))
      {
        G4double xa  = anActive->GetEnergy(a);
        theMerge->SetData(m, xa, anActive->GetXsec(a)+G4std::max(0., aPassive->GetXsec(xa)) );
        m++;
        a++;
        G4double xp = aPassive->GetEnergy(p);
        if( abs(abs(xp-xa)/xa)<0.001 )
        {
          p++;
        }
      } else {
        tmp = anActive; t=a;
        anActive = aPassive; a=p;
        aPassive = tmp; p=t;
      }
    }
    while (a!=anActive->GetVectorLength())
    {
      theMerge->SetData(m++, anActive->GetEnergy(a), anActive->GetXsec(a));
      a++;
    }
    while (p!=aPassive->GetVectorLength())
    {
      if(abs(theMerge->GetEnergy(G4std::max(0,m-1))-aPassive->GetEnergy(p))/aPassive->GetEnergy(p)>0.001)
        theMerge->SetData(m++, aPassive->GetEnergy(p), aPassive->GetXsec(p));
      p++;
    }
    delete theStore;
    theStore = theMerge;
  }

  G4ParticleChange * G4NeutronHPChannel::ApplyYourself(const G4Track & theTrack, G4int anIsotope)
  {
//    G4cout << "G4NeutronHPChannel::ApplyYourself+"<<niso<<G4endl;
    if(anIsotope != -1) return theFinalStates[anIsotope]->ApplyYourself(theTrack);
    G4double sum=0;
    G4int it=0;
    G4double * xsec = new G4double[niso];
    for (G4int i=0; i<niso; i++)
    {
      if(theFinalStates[i]->HasAnyData())
      {
        xsec[i] = theIsotopeWiseData[i].GetXsec(theTrack.GetKineticEnergy());
        sum += xsec[i];
      }
      else
      {
        xsec[i]=0;
      }
    } 
    if(sum == 0) 
    {
//      G4cout << "G4NeutronHPChannel::ApplyYourself theFinalState->Initialize+"<<G4endl;
//      G4cout << "G4NeutronHPChannel::ApplyYourself theFinalState->Initialize-"<<G4endl;
      it = niso*G4UniformRand();
    }
    else
    {
//      G4cout << "Are we still here? "<<sum<<G4endl;
//      G4cout << "TESTHP 23 NISO="<<niso<<G4endl;
      G4double random = G4UniformRand();
      G4double running=0;
//      G4cout << "G4NeutronHPChannel::ApplyYourself Done the sum"<<niso<<G4endl;
//      G4cout << "TESTHP 24 NISO="<<niso<<G4endl;
      for (G4int ix=0; ix<niso; ix++)
      {
        running += xsec[ix];
        if(random<=running/sum) 
        { 
          it = ix;
          goto OUT;
        }
      }
      OUT:
      if(it==niso) it--;
    }
    delete [] xsec;
    G4ParticleChange * theFinalState=NULL;
    while(theFinalState==NULL)
    {
//    G4cout << "TESTHP 24 it="<<it<<G4endl;
      theFinalState = theFinalStates[it]->ApplyYourself(theTrack);
    }
//    G4cout <<"THE IMPORTANT RETURN"<<G4endl;
    return theFinalState;
  }

