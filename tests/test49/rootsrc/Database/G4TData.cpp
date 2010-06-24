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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TData.h"


ClassImp(G4TData)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
void G4TData::CalcKinValues() // Calculates ProgMom/En & TargMass (in MeV)
{
  cout<<"in->G4TData::CalcKinValues: tg="<<fHeader.fTargetPDG<<", pr="
      <<fHeader.fProjectilePDG<<", Type="<<fHeader.fTypeVar<<", Val="<<fHeader.fTypeValue
      <<", Unit="<<fHeader.fTypeUnits<<endl;
  fTgMass         = gParticlesDAL->GetParticleMass(fHeader.fTargetPDG);
  Double_t prMass = gParticlesDAL->GetParticleMass(fHeader.fProjectilePDG);
  Double_t prM2   = prMass*prMass;
  Double_t prVal  = fHeader.fTypeValue;
  switch(fHeader.fTypeUnits)
  {
  case NoDim:
  case Degrees:
  case Radians:
    cout<<"Error>G4TData::CalcKinValues: Projectiles must have Emergy/Momentum Unit"<<endl;
    break;
  case keV:
    prVal/=1000.;
    break;
  case MeV:
    break;
  case GeV:
    prVal*=1000.;
    break;
  case TeV:
    prVal*=1000000.;
    break;
  }
  cout<<"Temporary>G4TData::CalcKinValues: EP = "<<prVal<<" MeV(/c), prM = "<<prMass<<endl;
  fPrMomentum=-1.;
  fPrEnergy=-1.;
  // Convert the value to total energy in MeV
  switch(fHeader.fTypeVar)
  {
  case Theta:
  case CosTheta:
  case LogTgHalfTheta:
  case p_L:
  case p_T:
  case x_F:
  case y_R:
    cout<<"Error>G4TData::CalcKinValues:Projectile must be Emergy/Momentum Variable"<<endl;
    break;
  case p_mom:
    fPrMomentum=prVal;
    fPrEnergy=std::sqrt(prVal*prVal+prM2);
    cout<<"Temporary>G4TData::CalcKinValues:(P): E="<<fPrEnergy<<", P="<<fPrMomentum<<endl;
    break;
  case E_tot:
    fPrEnergy=prVal;
    if(prVal<prMass)cout<<"Error>T?>G4TData::CalcKinValues:E="<<prVal<<"<M="<<prMass<<endl;
    fPrMomentum=std::sqrt(prVal*prVal-prM2);
    cout<<"Temporary>G4TData::CalcKinValues:(E): E="<<fPrEnergy<<", P="<<fPrMomentum<<endl;
    break;
  case T_kin:
    prVal+=prMass;
    fPrEnergy=prVal;
    fPrMomentum=std::sqrt(prVal*prVal-prM2);
    cout<<"Temporary>G4TData::CalcKinValues:(T): E="<<fPrEnergy<<", P="<<fPrMomentum<<endl;
    break;
  }  
  fCrossSection  = fHeader.fSigValue;
  switch(fHeader.fSigUnits)
  {
  case mBarn:
    break;
  case Barn:
    fCrossSection*=1000000.;
    break;
  case mkBarn:
    fCrossSection/=1000.;
    break;
  }
  cout<<"out->G4TData::CalcKinValues: prMom="<<fPrMomentum<<" MeV/c, prEn="<<fPrEnergy
      <<" MeV, tgMass="<<fTgMass<<" MeV, Sigma="<<fCrossSection<<" mb"<<endl;
}

//______________________________________________________________________________
void G4TData::Save()
{
  TString hfname = fDirectory + GetHeader() + ".root";
  //cout << "Saving the DataObject to file " << hfname << "..." << endl;
  // prepare the file
  TFile hfile(hfname, "recreate");
  if(hfile.IsZombie()) return;
  TObjArray hlist(0); // array to hold the list of objects
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    TString subHeader = fItems[i]->GetHeader();
    TTree* data = fItems[i]->GetData();
    hlist.Add(data);
  }
  TList* paramsList = new TList();
  TNamed* CrossSection = new TNamed("CrossSection",
                                    TString::Format("%g", fCrossSection).Data());
  TNamed* NumberOfEvents = new TNamed("NumberOfEvents",
                                      TString::Format("%d", fNumberOfEvents).Data());
  paramsList->Add(CrossSection);
  paramsList->Add(NumberOfEvents);
  hfile.WriteObject(paramsList,"Parameters","");

  hlist.Write();
  hfile.Close();

  TString name = this->GetHeader() + ".root";
  gCatalog->Insert(name);
}

//______________________________________________________________________________
void G4TData::Load(Int_t secondaryPDG)
{
  TString hfname = fDirectory + HeaderToString(fHeader) + ".root";
  TFile* hfile = new TFile(hfname.Data());
  if(hfile->IsZombie()) return;
  TString secName;
  if(secondaryPDG != 0) secName = gParticlesDAL->GetParticleName(secondaryPDG);
  TIter next(hfile->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next()))
  {
    TString name(key->GetName());
    if(name == "Parameters")
    {
      TList* paramsList = (TList*)hfile->Get(name.Data());
      fCrossSection  =  atof(((TNamed*)paramsList->At(0))->GetTitle());
      fNumberOfEvents =  atoi(((TNamed*)paramsList->At(1))->GetTitle());
      continue;
    }
    if(secondaryPDG != 0 && !name.BeginsWith(secName)) continue;
    G4TDataItem* item = this->AddItem(name);
    item->SetData( (TTree*) hfile->Get(name.Data()) );
  }
  fIsLoaded = true;
}

//______________________________________________________________________________
G4TDataItem* G4TData::AddItem(Int_t     SecondaryParticlePDG,
                              ArgEnum   CutVar,
                              UnitsEnum CutUnits,
                              Double_t  CutValue,
                              Double_t  CutDelta,
                              FunctEnum FunctionVar,
                              FunUnEnum FunctionUnits,
                              ErrorType FunctErrType,
                              ArgEnum   ArgumentVar,
                              UnitsEnum ArgumentUnits)
{
  G4TDataItem* item = new G4TDataItem(SecondaryParticlePDG,
                                      CutVar,
                                      CutUnits,
                                      CutValue,
                                      CutDelta,
                                      FunctionVar,
                                      FunctionUnits,
                                      FunctErrType,
                                      ArgumentVar,
                                      ArgumentUnits,
                                      fPrMomentum,
                                      fPrEnergy,
                                      fTgMass);
  fItems.push_back(item);
  return item;
}


//______________________________________________________________________________
G4TDataItem* G4TData::AddItem(TString const& headerStr)
{
  G4TDataItem* item = new G4TDataItem(headerStr, fPrMomentum, fPrEnergy, fTgMass);
  fItems.push_back(item);
  return item;
}

//______________________________________________________________________________
TString G4TData::HeaderToString(DataObjectHeader_t header) const
{
  // Be consistent with the enum's of G4TDataItem and this->StringToHeader()
  static const Int_t ntype=3;
  static const TString typ[ntype]={"Tkin" , "Pmom" , "Etot"};
  static const Int_t ntp=ntype-1;
  static const Int_t nunit=4;
  static const TString uni[nunit]={"MeV", "GeV", "TeV", "keV"};
  static const Int_t nun=nunit-1;
  if(header.fTypeVar > ntp || header.fTypeUnits > nun)
  {
    cout<<"Warning>G4TData::HeaderToString: Type="<<header.fTypeUnits<<" and/or Unit="
        <<header.fTypeVar<<" are beyon limits"<<endl;
    return "";
  }
  TString tp=typ[header.fTypeVar];
  TString un=uni[header.fTypeUnits];
  TString result = TString::Format("%s_%s_%g_%s_%s_%s_%s",
    gParticlesDAL->GetParticleName(header.fProjectilePDG).Data(), tp.Data(),
    header.fTypeValue, un.Data(), gParticlesDAL->GetParticleName(header.fTargetPDG).Data(),
    header.fModelName.Data(), header.fComment.Data());
  cout<<"G4TData::HeaderToString(out): Header ="<<result<<endl;
  return result;
}

//______________________________________________________________________________
DataObjectHeader_t G4TData::StringToHeader(TString headerStr) const
{
  // Be consistent with the enum's of G4TDataItem and this->HeaderToString()
  static const Int_t ntype=3;
  static const TString typ[ntype]={"Tkin" , "Pmom" , "Etot"};
  static const Int_t nunit=4;
  static const TString uni[nunit]={"MeV", "GeV", "TeV", "keV"};
  DataObjectHeader_t header;
  TString workStr(headerStr.ReplaceAll(".root",""));
  TObjArray* tokens = workStr.Tokenize("_");
  for(int i = 0; i < tokens->GetEntriesFast(); ++i)
  {
    TObjString *T = (TObjString*) (*tokens)[i];
    TString value = T->GetString();
    if     (!i)
    {
      header.fProjectilePDG  = gParticlesDAL->GetPDG(value.Data());
      cout<<"G4TData::StringToHeader: ProjPDG ="<<header.fProjectilePDG<<endl;
    }
    else if(i == 1)
    {
      TString tp=value.Data();
      if     (tp==typ[0]) header.fTypeVar=T_kin;
      else if(tp==typ[1]) header.fTypeVar=p_mom;
      else if(tp==typ[2]) header.fTypeVar=E_tot;
      else  cout<<"***Warning>G4TData::StringToHeader: ProjEnergy Type ="<<tp<<endl;
      cout<<"G4TData::StringToHeader: ProjEnergy Type ="<<header.fTypeVar<<endl;
    }
    else if(i == 2)
    {
      header.fTypeValue = atof(value.Data());
      cout<<"G4TData::StringToHeader: ProjEnergy Value ="<<header.fTypeValue<<endl;
    }
    else if(i == 3)
    {
      TString un=value.Data();
      if     (un==uni[0]) header.fTypeUnits=MeV;
      else if(un==uni[1]) header.fTypeUnits=GeV;
      else if(un==uni[2]) header.fTypeUnits=TeV;
      else if(un==uni[3]) header.fTypeUnits=keV;
      else  cout<<"***Warning>G4TData::StringToHeader: ProjEnergy Unit ="<<un<<endl;
      cout<<"G4TData::StringToHeader: ProjEnergy Units ="<<header.fTypeUnits<<endl;
    }
    else if(i == 4)
    {
      header.fTargetPDG = gParticlesDAL->GetPDG(value.Data());
      cout<<"G4TData::StringToHeader: Target PFG ="<<header.fTargetPDG<<endl;
    }
    else if(i == 5)
    {
      header.fModelName = value.Data();
      cout<<"G4TData::StringToHeader: Model Name ="<<header.fModelName<<endl;
    }
    else if(i == 6)
    {
      header.fComment = value.Data();
      cout<<"G4TData::StringToHeader: Comment ="<<header.fComment<<endl;
    }
  }
  header.fSigValue = 0.;
  header.fSigUnits = mBarn;
  return header;
}

//______________________________________________________________________________
void G4TData::BookHistograms(Int_t    hnbin,
                                Double_t hlxmin,
                                Double_t hlxmax,
                                Int_t    particleIdx,
                                Int_t    additionalIndex)
{
  // Keep Histograms in memory
  gROOT->cd();
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Double_t a   = item->GetCutValue();
    TString hname  = item->GetHistogramName(additionalIndex);
    TString htitle = TString::Format("T_{n} at %g deg (MeV)", a );
    cout<<"--->>G4TData::PrepareHistograms:Creating histogram "<<hname.Data()<<" at angle="
        <<a<<", nb="<<hnbin<<", xli="<<hlxmin<<", xla="<<hlxmax<<endl;
    TObject* existing = gROOT->FindObject(hname.Data());
    if(existing != 0) delete[] existing;
    TH1F* hist = new TH1F(hname.Data(), htitle.Data(), hnbin, hlxmin, hlxmax);
    item->SetHistogram(hist);
  }
}

//______________________________________________________________________________
vector<Double_t> G4TData::GetCutValues()
{
  vector<Double_t> resultVector;
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Double_t cutValue = item->GetCutValue();
    vector<Double_t>::iterator result;
    result = find( resultVector.begin(), resultVector.end(), cutValue );
    if(result == resultVector.end()) resultVector.push_back(cutValue);// IfNotFound->Insert
  }
  return resultVector;
}


//______________________________________________________________________________
vector<Int_t> G4TData::GetSecondaryPDGs()
{
  vector<Int_t> resultVector;
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Int_t pdg = item->GetSecondaryParticlePDG();
    vector<Int_t>::iterator result;
    result = find(resultVector.begin(), resultVector.end(), pdg);
    if(result == resultVector.end()) resultVector.push_back(pdg); // If not found, insert
  }
  return resultVector;
}

//______________________________________________________________________________
vector<Double_t> G4TData::GetCutValuesForSecondary(Int_t secondaryPDG)
{
  vector<Double_t> resultVector;
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Int_t pdg = item->GetSecondaryParticlePDG();
    if(pdg == secondaryPDG)
    {
      Double_t cutValue = item->GetCutValue();
      vector<Double_t>::iterator result;
      result = find( resultVector.begin(), resultVector.end(), cutValue );
      if(result == resultVector.end()) resultVector.push_back(cutValue);// notFound->Insert
    }
  }
  return resultVector;
}

//______________________________________________________________________________
vector<G4TDataItem*> G4TData::GetItemsForSecondary(Int_t secondaryPDG)
{
  cout<<"G4TData::GetItemsForSecondary:PDG="<<secondaryPDG<<",#ofIt="<<fItems.size()<<endl;
  vector<G4TDataItem*> resultVector;
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Int_t pdg = item->GetSecondaryParticlePDG();
    cout<<"G4TData::GetItemsForSecondary: PDG[="<<i<<"]="<<pdg<<endl;
    if(pdg == secondaryPDG)
    {
      vector<G4TDataItem*>::iterator result;
      result = find( resultVector.begin(), resultVector.end(), item );
      if(result == resultVector.end()) // If not found, insert
      {
        cout<<"G4TData::GetItemsForSecondary: Insert NewPDG="<<pdg<<endl;
        resultVector.push_back(item);
      }
    }
  }
  return resultVector;
}

//______________________________________________________________________________
vector<G4TDataItem*> G4TData::GetItemsForSecondary(vector<Int_t> secondaries)
{
  vector<G4TDataItem*> resultVector;
  for(UInt_t i = 0; i < secondaries.size(); ++i)
  {
    vector<G4TDataItem*> items = GetItemsForSecondary(secondaries[i]);
    for(UInt_t j = 0; j < items.size(); ++j) resultVector.push_back(items[j]);
  }
  return resultVector;
}

//______________________________________________________________________________
G4TDataItem* G4TData::GetItem(Int_t secondaryPDG, Double_t cutValue)
{
  G4TDataItem* result = 0;
  for(UInt_t i = 0; i < fItems.size(); ++i)
  {
    G4TDataItem* item = fItems[i];
    Int_t pdg = item->GetSecondaryParticlePDG();
    Double_t cut = item->GetCutValue();
    if(pdg == secondaryPDG && cut == cutValue)
    {
      result = item;
      break;
    }
  }
  return result;
}

//______________________________________________________________________________
std::pair<Double_t, Double_t> G4TData::GetLimits(Int_t secondaryIdx, Int_t padsPerRow)
{
  cout<<"G4TData::GetLimits: ID="<<secondaryIdx<<", #0fPads="<<padsPerRow<<endl;
  Double_t rMin=0.;
  Double_t rMax=0.;
  vector<Int_t> fragments = GetSecondaryPDGs();
  Bool_t initialized = false;
  // Auto scale feature
  if(secondaryIdx % padsPerRow == 0)
  {
    UInt_t max = secondaryIdx + padsPerRow;
    cout<<"G4TData::GetLimits: max="<<max<<endl;
    if(max > fragments.size()) max = fragments.size();
    for(UInt_t k = secondaryIdx; k < max; ++k)
    {
      Int_t pdgCode = fragments[k];
      vector<G4TDataItem*> citems = GetItemsForSecondary(pdgCode);
      cout<<"G4TData::GetLimits:k="<<k<<", PDG="<<pdgCode<<", #ofIt="<<citems.size()<<endl;
      for(UInt_t l = 0; l < citems.size(); ++l)
      {
        cout<<"G4TData::GetLimits: k="<<k<<", l="<<l<<",pointerToIt="<<citems[l]<<endl;
        std::pair<Double_t, Double_t> limits = citems[l]->GetLimits();
        Double_t cMin=limits.first;
        Double_t cMax=limits.second;
        cout<<"G4TData::GetLimits: l="<<l<<", cMin="<<cMin<<", cMax="<<cMax<<endl;
        if(!initialized)
        {
          rMax = cMax;
          rMin = cMin;
          initialized = true;
        }
        else
        {
          if(cMax > rMax ) rMax = cMax;
          if(cMin < rMin ) rMin = cMin;
        }
        cout<<"G4TData::GetLimits: l="<<l<<", rMin="<<rMin<<", rMax="<<rMax<<endl;
      }
    }
  }
  cout<<">>>G4TData::GetLimits(out): min = "<<rMin<<", max = "<< rMax<<endl;
  return std::make_pair(rMin,rMax);
}

//______________________________________________________________________________
std::pair<Double_t, Double_t> G4TData::GetLimits()
{
  Double_t rMin=0.;
  Double_t rMax=0.;
  vector<Int_t> fragments = GetSecondaryPDGs();
  Bool_t initialized = false;
  // Auto Y scale feature
  for(UInt_t k = 0; k < fragments.size(); ++k)
  {
    Int_t pdgCode = fragments[k];
    vector<G4TDataItem*> citems = this->GetItemsForSecondary(pdgCode);
    for(UInt_t l = 0; l < citems.size(); ++l)
    {
      std::pair<Double_t, Double_t> limits = citems[l]->GetLimits();
      Double_t cMin=limits.first;
      Double_t cMax=limits.second;
      if(!initialized)
      {
        rMax = cMax;
        rMin = cMin;
        initialized = true;
      }
      else
      {
        if(cMax > rMax ) rMax = cMax;
        if(cMin < rMin ) rMin = cMin;
      }
    }
  }
  cout<<">>>G4TData::GetLimits(): minY = "<<rMin<<", maxY = "<< rMax<<endl;
  return std::make_pair(rMin,rMax);
}

//______________________________________________________________________________
Double_t G4TData::GetMaxT(Int_t secondaryIdx, Int_t padsPerRow)
{
  cout<<"G4TData::GetMaxT: ID="<<secondaryIdx<<", #0fPads="<<padsPerRow<<endl;
  Double_t maxT=0.;
  vector<Int_t> fragments = GetSecondaryPDGs();
  Bool_t initialized = false;
  // Auto T_kin scale feature
  if(secondaryIdx % padsPerRow == 0)
  {
    UInt_t max = secondaryIdx + padsPerRow;
    cout<<"G4TData::GetMaxT: max="<<max<<endl;
    if(max > fragments.size()) max = fragments.size();
    for(UInt_t k = secondaryIdx; k < max; ++k)
    {
      Int_t pdgCode = fragments[k];
      vector<G4TDataItem*> citems = GetItemsForSecondary(pdgCode);
      cout<<"G4TData::GetMaxT:k="<<k<<", PDG="<<pdgCode<<", #ofIt="<<citems.size()<<endl;
      for(UInt_t l = 0; l < citems.size(); ++l)
      {
        cout<<"G4TData::GetMaxT: k="<<k<<", l="<<l<<",pointerToIt="<<citems[l]<<endl;
        Double_t cMax = citems[l]->GetMaxT();
        cout<<"G4TData::GetMaxT: l="<<l<<", cMax="<<cMax<<endl;
        if(!initialized)
        {
          maxT = cMax;
          initialized = true;
        }
        else if(cMax > maxT ) maxT = cMax;
        cout<<"G4TData::GetMaxT: l="<<l<<", maxT="<<maxT<<endl;
      }
    }
  }
  cout<<">>>G4TData::GetMaxT(out): maxT = "<< maxT<<endl;
  return maxT*1.03;
}

//______________________________________________________________________________
Double_t G4TData::GetMaxT()
{
  Double_t maxT=0.;
  vector<Int_t> fragments = GetSecondaryPDGs();
  Bool_t initialized = false;
  // Auto T_kin scale feature
  for(UInt_t k = 0; k < fragments.size(); ++k)
  {
    Int_t pdgCode = fragments[k];
    vector<G4TDataItem*> citems = this->GetItemsForSecondary(pdgCode);
    for(UInt_t l = 0; l < citems.size(); ++l)
    {
      Double_t cMax = citems[l]->GetMaxT();
      if(!initialized)
      {
        maxT = cMax;
        initialized = true;
      }
      else if(cMax > maxT ) maxT = cMax;
    }
  }
  cout<<">>>G4TData::MaxT(): maxT = "<< maxT<<endl;
  return maxT;
}
