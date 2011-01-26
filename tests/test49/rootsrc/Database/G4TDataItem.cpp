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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TDataItem.h"

ClassImp(G4TDataItem)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
G4TDataItem::G4TDataItem(TString const& headerStr, Double_t  PrMomentum,
                         Double_t  PrEnergy, Double_t  TgMass) // XS=0.
{
  fHeader     = StringToHeader(headerStr);
  fHistogram  = 0;
  fData       = 0;
  fPrMomentum = PrMomentum;
  fPrEnergy   = PrEnergy;
  fTgMass     = TgMass;
  cout<<"G4TDataItem::Constr(H): PDG="<<fHeader.fSecondaryParticlePDG<<", CT="
      <<fHeader.fCutVar<<", CU="<<fHeader.fCutUnits<<", CV="<<fHeader.fCutValue<<", CD="
      <<fHeader.fCutDelta<<", FT="<<fHeader.fFunctionVar<<", FU="<<fHeader.fFunctionUnits
      <<", FE="<<fHeader.fFunctErrorType<<", AT="<<fHeader.fArgumentVar<<", AU="
      <<fHeader.fArgumentUnits<<", pP="<<fPrMomentum<<", pE="<<fPrEnergy<<", tM="<<fTgMass
      <<endl;
}

//______________________________________________________________________________
G4TDataItem::G4TDataItem(Int_t     SecondaryParticlePDG,
                         ArgEnum   CutVar,
                         UnitsEnum CutUnits,
                         Double_t  CutValue,
                         Double_t  CutDelta,
                         FunctEnum FunctionVar,
                         FunUnEnum FunctionUnits,
                         ErrorType FunctErrorType,
                         ArgEnum   ArgumentVar,
                         UnitsEnum ArgumentUnits,
                         Double_t  PrMomentum,
                         Double_t  PrEnergy,
                         Double_t  TgMass)
{
  fHeader.fSecondaryParticlePDG = SecondaryParticlePDG;
  fHeader.fCutVar               = CutVar;
  fHeader.fCutUnits             = CutUnits;
  fHeader.fCutValue             = CutValue;
  fHeader.fCutDelta             = CutDelta;
  fHeader.fFunctionVar          = FunctionVar;
  fHeader.fFunctionUnits        = FunctionUnits;
  fHeader.fFunctErrorType       = FunctErrorType;
  fHeader.fArgumentVar          = ArgumentVar;
  fHeader.fArgumentUnits        = ArgumentUnits;
  fHistogram                    = 0;
  fData                         = 0;
  fPrMomentum                   = PrMomentum;
  fPrEnergy                     = PrEnergy;
  fTgMass                       = TgMass;
  cout<<"G4TDataItem::Constr(D): PDG="<<fHeader.fSecondaryParticlePDG<<", CT="
      <<fHeader.fCutVar<<", CU="<<fHeader.fCutUnits<<", CV="<<fHeader.fCutValue<<", CD="
      <<fHeader.fCutDelta<<", FT="<<fHeader.fFunctionVar<<", FU="<<fHeader.fFunctionUnits
      <<", FE="<<fHeader.fFunctErrorType<<", AT="<<fHeader.fArgumentVar<<", AU="
      <<fHeader.fArgumentUnits<<", pP="<<fPrMomentum<<", pE="<<fPrEnergy<<", tM="<<fTgMass
      <<endl;
}

//______________________________________________________________________________
void G4TDataItem::LoadFromASCII(const TString& filename, Bool_t load3)
{
  gROOT->cd();
  TString tname = GetHeader();
  TString vname = "t:s";                 // name set of the read out values (default:T,Sig)
  fData = new TTree(tname.Data(), tname.Data());
  Float_t dval = 0;
  Float_t sval = 0;
  Float_t rval = 0; // The branch "r" is not used later on (@@ how to get rid of it?)
  //
  TBranch* branchS = 0;
  TBranch* branchD = 0;
  Int_t   rEntries = 0;
  // @@ Why switch does not work on some compilers ?
  //switch(fHeader.fFunctErrorType)
  //{
  //case NoError:
  if(fHeader.fFunctErrorType == NoError)
  {
    fData->ReadFile(filename.Data(),vname.Data());
    // wait for tree
    while(fData == NULL) gSystem->Sleep(10);
    // @@ OPTIMIZE (need a better way to do this)
    // in kumac: ve/cop $SIGMA(s[j][a]/20.) d[j][a]
    branchS = fData->GetBranch("s");
    branchD = fData->Branch("d", &dval, "d");
    TBranch* branchS = fData->GetBranch("s");
    TBranch* branchD = fData->Branch("d", &dval, "d");
    branchS->SetAddress(&sval);
    rEntries = fData->GetEntries();
    //Int_t rEntries = fData->GetEntries();
    for(Int_t i = 0; i < rEntries; ++i)
    {
      branchS->GetEntry(i);
      dval = sval / 20;                  // Fake 5% errors
      branchD->Fill();
    }
  }
  else if(fHeader.fFunctErrorType == Absolute)
  {
		//  break;
  //case Absolute:
    vname = "t:s:d";
    fData->ReadFile(filename.Data(),vname.Data());
  }
  else if(fHeader.fFunctErrorType == Absolute)
  {
  //  break;
  //case Relative:
    vname = "t:s:r";
    fData->ReadFile(filename.Data(),vname.Data());
    // wait for tree
    while(fData == NULL) gSystem->Sleep(10);
    // @@ OPTIMIZE (need a better way to do this)
    // in kumac: ve/cop $SIGMA(r1[a]*s1[a]/100.) d1[a]
    TBranch* branchR = fData->GetBranch("r");
    branchS = fData->GetBranch("s");
    branchD = fData->Branch("d", &dval, "d");
    branchR->SetAddress(&rval);
    branchS->SetAddress(&sval);
    rEntries = fData->GetEntries();
    for(Int_t i = 0; i < rEntries; ++i)
    {
      branchR->GetEntry(i);
      branchS->GetEntry(i);
      dval = rval * sval;
      branchD->Fill();
    }
  }
  else if(fHeader.fFunctErrorType == InPerCent)
  {
  //  break;
  //case InPerCent:
    vname = "t:s:p";
    fData->ReadFile(filename.Data(),vname.Data());
    // wait for tree
    while(fData == NULL) gSystem->Sleep(10);
    // @@ OPTIMIZE (need a better way to do this)
    // in kumac: ve/cop $SIGMA(r1[a]*s1[a]/100.) d1[a]
    TBranch* branchP = fData->GetBranch("p");
    branchS = fData->GetBranch("s");
    branchD = fData->Branch("d", &dval, "d");
    branchP->SetAddress(&rval);
    branchS->SetAddress(&sval);
    rEntries = fData->GetEntries();
    for(Int_t i = 0; i < rEntries; ++i)
    {
      branchP->GetEntry(i);
      branchS->GetEntry(i);
      dval = rval * sval / 100;
      branchD->Fill();
    }
  }
  else
  {
  //  break;
  //default:
    cout<<"*Warning*G4TDataItem::LoadFromASCII:errorType= "<<fHeader.fFunctErrorType<<endl;
  }
}

//______________________________________________________________________________
TString G4TDataItem::HeaderToString(DataItemObjectHeader_t header) const
{
  cout<<"G4TDataItem::HeaderToString(in): PDG="<<fHeader.fSecondaryParticlePDG<<", CT="
      <<fHeader.fCutVar<<", CU="<<fHeader.fCutUnits<<", CV="<<fHeader.fCutValue<<", CD="
      <<fHeader.fCutDelta<<", FT="<<fHeader.fFunctionVar<<", FU="<<fHeader.fFunctionUnits
      <<", FE="<<fHeader.fFunctErrorType<<", AT="<<fHeader.fArgumentVar<<", AU="
      <<fHeader.fArgumentUnits<<", pP="<<fPrMomentum<<", pE="<<fPrEnergy<<", tM="<<fTgMass
      <<endl;
  TString result = TString::Format("%s_%d_%d_%g_%g_%d_%d_%d_%d_%d",
    gParticlesDAL->GetParticleName(fHeader.fSecondaryParticlePDG).Data(),
    fHeader.fCutVar,
    fHeader.fCutUnits,
    fHeader.fCutValue,
    fHeader.fCutDelta,
    fHeader.fFunctionVar,
    fHeader.fFunctionUnits,
    fHeader.fFunctErrorType,
    fHeader.fArgumentVar,
    fHeader.fArgumentUnits);
  cout<<"G4TDataItem::HeaderToString(out): "<<result<<endl;
  return result;
}

//______________________________________________________________________________
DataItemObjectHeader_t G4TDataItem::StringToHeader(TString headerStr) const
{
  cout<<"G4TDataItem::StringToHeader(in): "<<headerStr<<endl;
  DataItemObjectHeader_t header;
  TObjArray* tokens = headerStr.Tokenize("_");
  for(Int_t i = 0; i< tokens->GetEntriesFast(); ++i)
  {
    TObjString* T = (TObjString*) (*tokens)[i];
    TString value = T->GetString();
    if(i == 0)
    {
      cout<<"G4TDataItem::StringToHeader(PDG): i="<<i<<", v="<<value.Data()<<endl;
      header.fSecondaryParticlePDG = gParticlesDAL->GetPDG(value.Data());
    }
    else if(i == 1)
    {
      cout<<"G4TDataItem::StringToHeader(CV): i="<<i<<", v="<<value.Data()<<endl;
      header.fCutVar               = (ArgEnum)atoi(value.Data());
    }
    else if(i == 2)
    {
      cout<<"G4TDataItem::StringToHeader(CU): i="<<i<<", v="<<value.Data()<<endl;
      header.fCutUnits             = (UnitsEnum)atoi(value.Data());
    }
    else if(i == 3)
    {
      cout<<"G4TDataItem::StringToHeader(CV): i="<<i<<", v="<<value.Data()<<endl;
      header.fCutValue             = atof(value.Data());
    }
    else if(i == 4)
    {
      cout<<"G4TDataItem::StringToHeader(CD): i="<<i<<", v="<<value.Data()<<endl;
      header.fCutDelta             = atof(value.Data());
    }
    else if(i == 5)
    {
      cout<<"G4TDataItem::StringToHeader(FT): i="<<i<<", v="<<value.Data()<<endl;
      header.fFunctionVar          = (FunctEnum)atoi(value.Data());
    }
    else if(i == 6)
    {
      cout<<"G4TDataItem::StringToHeader(FU): i="<<i<<", v="<<value.Data()<<endl;
      header.fFunctionUnits        = (FunUnEnum)atoi(value.Data());
    }
    else if(i == 7)
    {
      cout<<"G4TDataItem::StringToHeader(FE): i="<<i<<", v="<<value.Data()<<endl;
      header.fFunctErrorType       = (ErrorType)atoi(value.Data());
    }
    else if(i == 8)
    {
      cout<<"G4TDataItem::StringToHeader(AT): i="<<i<<", v="<<value.Data()<<endl;
      header.fArgumentVar          = (ArgEnum)atoi(value.Data());
    }
    else if(i == 9)
    {
      cout<<"G4TDataItem::StringToHeader(AU): i="<<i<<", v="<<value.Data()<<endl;
      header.fArgumentUnits        = (UnitsEnum)atoi(value.Data());
    }
  }
  cout<<"G4TDataItem::StringToHeader(out): PDG="<<header.fSecondaryParticlePDG<<", CT="
      <<header.fCutVar<<", CU="<<header.fCutUnits<<", CV="<<header.fCutValue<<", CD="
      <<header.fCutDelta<<", FT="<<header.fFunctionVar<<", FU="<<header.fFunctionUnits
      <<", FE="<<header.fFunctErrorType<<", AT="<<header.fArgumentVar<<", AU="
      <<header.fArgumentUnits<<endl;
  return header;
}

//______________________________________________________________________________
//TString G4TDataItem::MakeInvariantFunction()
//{
//  Double_t mass = gParticlesDAL->GetParticleMass(fHeader.fSecondaryParticlePDG);
//  TString P = TString::Format("(sqrt(t * (2. * %g + t)))", mass);
//  TString DFormula = TString::Format("(s/%s) : t : (d/%s)", P.Data(), P.Data());
//  return DFormula;
//}

//____________________________________________________ @@ make similar for 1/4piA
TString G4TDataItem::GetInvariantFunction()  // Calculates momentum(p) & devides s & d by p
{
  Double_t  aW = GetTargMassMeV()/931.494043;// A (atomic Weight) in UAMU
  FunctEnum funVar = GetFunctionVar();       // Type of the function
  FunUnEnum funUnit= GetFunctionUnits();     // Units of the function
  ArgEnum   argVar = GetArgumentVar();       // Type of the Argument
  UnitsEnum argUnit= GetArgumentUnits();     // Units of the Arguent
  Double_t mass = gParticlesDAL->GetParticleMass(fHeader.fSecondaryParticlePDG); // in MeV
  Double_t mass2= mass * mass;
  cout<<"G4TDataItem::GetLimits(in): A_t="<<aW<<", FT="<<funVar<<", FU="<<funUnit<<", AT="
      <<argVar<<", AU="<<argUnit<<", secMass="<<mass<<endl;
  // @@ T(=t) should be always in MeV (!) to correspond to the global units (!) @@ SIGMA(?)
  TString P = "t"; // Momentum for normalization of the invariant function
  TString E = "t"; // Kinetic energy as the second column
  if(argVar==T_kin)
  {
    if     (argUnit==MeV)
    {
      P = TString::Format("(sqrt(t * (2. * %g + t)))", mass);
      cout<<"G4TDataItem::GetLimits(in): T_kin MeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==GeV)
    {
      P = TString::Format("(sqrt(t*1000.*(2. * %g + t*1000.)))", mass);
      E = "(t*1000.)";
      cout<<"G4TDataItem::GetLimits(in): T_kin GeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==TeV)
    {
      P = TString::Format("(sqrt(t*1000000.*(2. * %g + t*1000000.)))", mass);
      E = "(t*1000000.)";
      cout<<"G4TDataItem::GetLimits(in): T_kin TeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==keV)
    {
      P = TString::Format("(sqrt(t/1000.*(2. * %g + t/1000.)))", mass);
      E = "(t/1000.)";
      cout<<"G4TDataItem::GetLimits(in): T_kin keV P="<<P<<", E="<<E<<endl;
    }
    else cout<<"G4TDataItem::GetInvariantFunction: T_kin must be in kMGTeV"<<endl;
  }
  else if(argVar==E_tot)
  {
    if     (argUnit==MeV)
    {
      P = TString::Format("(sqrt(t ** 2 - %g))", mass2);
      E = TString::Format("(t - %g)", mass);
      cout<<"G4TDataItem::GetLimits(in): E_tot MeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==GeV)
    {
      P = TString::Format("(sqrt((t*1000.) ** 2. - %g))", mass2);
      E = TString::Format("(t*1000. - %g)", mass);
      cout<<"G4TDataItem::GetLimits(in): E_tot GeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==TeV)
    {
      P = TString::Format("(sqrt((t*1000000.)**2 - %g))", mass2);
      E = TString::Format("(t*1000000. - %g)", mass);
      cout<<"G4TDataItem::GetLimits(in): E_tot TeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==keV)
    {
      P = TString::Format("(sqrt(t/1000.)**2 - %g))", mass2);
      E = TString::Format("(t/1000. - %g)", mass);
      cout<<"G4TDataItem::GetLimits(in): E_tot keV P="<<P<<", E="<<E<<endl;
    }
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: E_tot must be in kMGTeV"<<endl;
  }
  else if(argVar==p_mom)
  {
    if     (argUnit==MeV)
    {
      E = TString::Format("(sqrt(t ** 2 + %g) - %g)", mass2, mass);
      cout<<"G4TDataItem::GetLimits(in): p_mom MeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==GeV)
    {
      E = TString::Format("(sqrt((t*1000.) ** 2. + %g) - %g)", mass2, mass);
      P = "(t*1000.)";
      cout<<"G4TDataItem::GetLimits(in): p_mom GeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==TeV)
    {
      E = TString::Format("(sqrt((t*1000000.)**2 + %g) - %g)", mass2, mass);
      P = "(t*1000000.)";
      cout<<"G4TDataItem::GetLimits(in): p_mom TeV P="<<P<<", E="<<E<<endl;
    }
    else if(argUnit==keV)
    {
      E = TString::Format("(sqrt(t/1000.)**2 + %g) - %g)", mass2, mass);
      P = "(t/1000.)";
      cout<<"G4TDataItem::GetLimits(in): p_mom keV P="<<P<<", E="<<E<<endl;
    }
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: p_mom must be in kMGTeV"<<endl;
  }
  else cout<<"Warning>G4TDataItem::GetInvariantFunction:EP="<<argVar<<" (use E,T,p)"<<endl;
  TString DFormula = "";
  if     (funVar==dN_dEdO) // @@ not multiplied by the sigma yet (Improve!)
  {
    if     (funUnit==_MeV_sr)
      DFormula = TString::Format("(s/%s) : %s : (d/%s)", P.Data(), E.Data(), P.Data());
    else if(funUnit==_GeV_sr)
      DFormula = TString::Format("(s*.001/%s) : %s : (d*.001/%s)",
                                 P.Data(), E.Data(), P.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dN_dEdO"<<endl;
  }
  else if(funVar==dS_dEdO)
  {
    if     (funUnit==mb_MeV_sr)
      DFormula = TString::Format("(s/%s) : %s : (d/%s)", P.Data(), E.Data(), P.Data());
    else if(funUnit==mb_GeV_sr)
      DFormula = TString::Format("(s/1000./%s) : %s : (d/1000./%s)",
                                 P.Data(), E.Data(), P.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dS_dEdO"<<endl;
  }
  else if(funVar==dN_dpdO) // @@ not multiplied by the sigma
  {
    if     (funUnit==_MeV_sr)
      DFormula = TString::Format("(s*%s/%s/%s) : %s : (d*%s/%s/%s)", E.Data(), P.Data(),
                                 P.Data(), E.Data(), E.Data(), P.Data(), P.Data());
    else if(funUnit==_GeV_sr)
      DFormula = TString::Format("(s*%s/1000./%s/%s) : %s : (d*%s/1000/%s/%s)", E.Data(),
                                 P.Data(),P.Data(), E.Data(), E.Data(),P.Data(),P.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dN_dpdO"<<endl;
  }
  else if(funVar==dS_dpdO)
  {
    if     (funUnit==mb_MeV_sr)
      DFormula = TString::Format("(s*%s/%s/%s) : %s : (d*%s/%s/%s)", E.Data(), P.Data(),
                                 P.Data(), E.Data(), E.Data(), P.Data(), P.Data());
    else if(funUnit==mb_GeV_sr)
      DFormula = TString::Format("(s*%s/1000./%s/%s) : %s : (d*%s/1000/%s/%s)", E.Data(),
                                 P.Data(),P.Data(), E.Data(), E.Data(),P.Data(),P.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dN_dpdO"<<endl;
  }
  else if(funVar==dN_pdEdO) // @@ not multiplied by the sigma
  {
    if     (funUnit==_MeV2_sr)
      DFormula = TString::Format("s : %s : d", E.Data());
    else if(funUnit==_GeV2_sr)
      DFormula = TString::Format("(s/1000000) : %s : (d/1000000)", E.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dN_pdEdO"<<endl;
  }
  else if(funVar==dS_pdEdO)
  {
    if     (funUnit==mb_MeV2_sr)
      DFormula = TString::Format("s : %s : d", E.Data());
    else if(funUnit==mb_GeV2_sr)
      DFormula = TString::Format("(s/1000000) : %s : (d/1000000)", E.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = dS_pdEdO"<<endl;
  }
  else if(funVar==EdN_d3p) // @@ not multiplied by sigma (same as dN_pdEdO, @@ can be ||)
  {
    if     (funUnit==_MeV2_sr)
      DFormula = TString::Format("s : %s : d", E.Data());
    else if(funUnit==_GeV2_sr)
      DFormula = TString::Format("(s/1000000) : %s : (d/1000000)", E.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = EdN_d3p"<<endl;
  }
  else if(funVar==EdS_d3p) // (same as dS_pdEdO, @@ can be ||)
  {
    if     (funUnit==mb_MeV2_sr)
      DFormula = TString::Format("s : %s : d", E.Data());
    else if(funUnit==mb_GeV2_sr)
      DFormula = TString::Format("(s/1000000) : %s : (d/1000000)", E.Data());
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = EdS_d3p"<<endl;
  }
  else if(funVar==EdS_Ad3p) // @@ not divided by the sigma & mom (same as dS_pdEdO)
  {
    cout<<"G4TDataItem::GetLimits:*FT=EdS_Ad3p* FU="<<funUnit<<", A="<<aW<<", P="<<P.Data()
        <<", E="<<E.Data()<<endl;
    if     (funUnit==mb_MeV2_sr)
      DFormula = TString::Format("(s*%g) : %s : (d*%g)",
                                 aW, E.Data(), aW);
    else if(funUnit==mb_GeV2_sr)
      DFormula = TString::Format("(s*%g/1000000) : %s : (d*%g/1000000)",
                                 aW, E.Data(), aW);
    else cout<<"Warning>G4TDataItem::GetInvariantFunction: Wrong units="<<funUnit
             <<" for Function type = EdS_Ad3p"<<endl;
  }
  else cout<<"***Warning>G4TDataItem::GetInvariantFunction: FunType="<<funVar
           <<" isn't yet supported"<<endl; // @@ the growing point for other functions
  cout<<"G4TDataItem::GetInvariantFunction(out): InvF="<<DFormula<<endl;
  return DFormula;
}

//________________________________________________________________________________________
Double_t G4TDataItem::GetMaxT() // Define the maximum kinetic energy
{
  Double_t maxT=0.;
  cout<<"G4TDataItem::GetMaxT(in): fData="<<fData<<endl;
  Bool_t   initialized = false;
  if(fData != 0)
  {
    fData->Draw(GetInvariantFunction().Data(),"","goff"); // Draw M.K.?
    cout<<"G4TDataItem::GetMaxT: Invariant function is Drown "<<endl;
    Long64_t nb = fData->GetSelectedRows();
    Double_t* p = fData->GetV2();           // pinter to the first column value
    cout<<"G4TDataItem::MaxT: #ofRows = "<<nb<<", firstInColumn2 = "<<*p<<endl;
    for(Int_t i = 0; i < nb; ++i)
    {
      if(!initialized) // the first point initializes the min & max values
      {
        maxT = *p;
        initialized = true;
      }
      else if(*p > maxT) maxT = *p;
      p++;
    }
  }
  else cout<<"-Warning->G4TDataItem::GetMaxT: fData=0, Data are not defined"<<endl;
  cout<<"G4TDataItem::GetMaxT(out): maxT = "<<maxT<<endl;
  return maxT;
}

//________________________________________________________________________________________
std::pair<Double_t, Double_t> G4TDataItem::GetLimits() // Define the range of the function
{
  Double_t rMin=0.;
  Double_t rMax=0.;
  cout<<"G4TDataItem::GetLimits(in): fData="<<fData<<endl;
  Bool_t   initialized = false;
  if(fData != 0)
  {
    fData->Draw(GetInvariantFunction().Data(),"","goff"); // Draw M.K.?
    cout<<"G4TDataItem::GetLimits: Invariant function is Drown "<<endl;
    Long64_t nb = fData->GetSelectedRows();
    Double_t* p = fData->GetV1();           // pinter to the first column value
    cout<<"G4TDataItem::GetLimits: #ofRows = "<<nb<<", vColumn1 = "<<*p<<endl;
    for(Int_t i = 0; i < nb; ++i)
    {
      if(!initialized) // the first point initializes the min & max values
      {
        rMin = *p;
        rMax = *p;
        initialized = true;
      }
      else
      {
        if(*p < rMin) rMin = *p;
        if(*p > rMax) rMax = *p;
      }
      p++;
    }
  }
  else cout<<"-Warning->G4TDataItem::GetLimits: fData=0, Data are not defined"<<endl;
  cout<<"G4TDataItem::GetLimits(out): min = "<<rMin<<", max = "<<rMax<<endl;
  return std::make_pair(rMin, rMax);
}

//______________________________________________________________________________
Double_t G4TDataItem::GetAngleInRadians() const
{
  // Set Angle
  Double_t  a      = GetCutValue();
  UnitsEnum cutUnit= GetCutUnits();
  ArgEnum   cutVar = GetCutVar();
  Double_t arad = 0;
  if      (cutVar == Theta)
  {
    if     (cutUnit == Radians) arad = a;
    else if(GetCutUnits() == Degrees ) arad = 3.14159265 * a / 180.;
    else cout<<"Warning>G4TDataItem::GetAngleInRadians:Angle NOT in degrees/radians"<<endl;
  }
  else if (cutVar == CosTheta)
  {
    if   (cutUnit == NoDim) arad = std::acos(a);
    else cout<<"Warning>G4TDataItem::GetAngleInRadians:Cos(Theta) should be NoDim"<<endl;
  }
  else if (cutVar == LogTgHalfTheta)
  {
    if   (cutUnit == NoDim) arad = 2*std::atan(std::exp(a));
    else cout<<"Warning>G4TDataItem::GetAngleInRadians: eta should be NoDim"<<endl;
  }
  else cout<<"Warning>G4TDataItem::GetAngleInRadians: *** Not angle vadiable ***"<<endl;
  return arad;
}

//______________________________________________________________________________
TString G4TDataItem::GetHistogramName(Int_t anIndex ) const
{
  if(anIndex == -1) return TString::Format("Histogram_%s", HeaderToString(fHeader).Data());
  return   TString::Format("Histogram_%s_%d", HeaderToString(fHeader).Data(), anIndex);
}
