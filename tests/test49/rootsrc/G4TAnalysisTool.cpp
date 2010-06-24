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

#include "G4TAnalysisTool.h"

G4TAnalysisTool *gAnalysisTool = new G4TAnalysisTool();

ClassImp(G4TAnalysisTool)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
Int_t G4TAnalysisTool::Run(TString const& publicationFile, Int_t secondaryPDGorIdx,
                           TString const& printQue )
{
  cout<<"G4TAnalysisTool::Run: Loading Libraries..."<<endl;
  G4TSimHelper::LoadLibraries();
  gBenchmark->Start("Overall Benchmark");
  gROOT->Reset();
  // parse the header and load the publication
  fPublication = gTestDB->LoadData(publicationFile);
  // prepare arguments
  fProjectilePDG  = fPublication->fHeader.fProjectilePDG;
  fTargetPDG   = fPublication->fHeader.fTargetPDG;
  fModelName  = fPublication->fHeader.fModelName;
  fKineticEnergy  = fPublication->fHeader.fTypeValue;
  // get nn and np
  fNN = gParticlesDAL->GetN(fTargetPDG);
  fNP = gParticlesDAL->GetZ(fTargetPDG);
  // initialize and execute
  Initialize();
  InternalExecute(secondaryPDGorIdx, fNP, fNN, fKineticEnergy, printQue, 2, fNP);

  gBenchmark->Show("Overall Benchmark");
  return 0;
}

//______________________________________________________________________________
void G4TAnalysisTool::InternalExecute(Int_t secondaryPDGIndx, Int_t np, Int_t nn,
                                      Double_t en, const TString& pq, Int_t nzone,
                                      Int_t nvex, const TString& dir)
{
  if(fPublication == 0 || !fPublication->IsLoaded())
  {
    cout<<"G4TAnalysisTool::InternalExecute: Publication isn't defined! Aborting..."<<endl;
    return;
  }
  // Redefine xmax and prepare general Hisograms
  fHxmax = fPublication->GetMaxT() * 1.03; // Step a bit to the right from the max point
  fHlxmax = std::log(fHxmax);
  cout<<"G4TAnalysisTool::InternalExecute:Cuts&Hists, xi="<<fHlxmin<<",xa="<<fHlxmax<<endl;
  PrepareHistograms(fHnbin, fHlxmin, fHlxmax); // Taken out of the G4TTool Constructor
  SetSecondaryToAnalyze(secondaryPDGIndx, np);
  // Prepare Models
  ArgEnum arg      = fPublication->fHeader.fTypeVar;
  UnitsEnum units  = fPublication->fHeader.fTypeUnits;
  Double_t value   = fPublication->fHeader.fTypeValue;
  SigmaEnum sigun  = fPublication->fHeader.fSigUnits;
  Double_t sigval  = fPublication->fHeader.fSigValue;
  
  fSimulations.clear();
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "lhep", "No", arg,
                         value, units, sigval, sigun,  1));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "chips", "No", arg,
                         value, units, sigval, sigun,  2));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "preco", "No", arg,
                         value, units, sigval, sigun,  6));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "binary", "No",
                         arg, value, units, sigval, sigun,  8));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "bertini", "No",
                         arg, value, units, sigval, sigun,  9));
  
  vector<G4TDataItem*> pubItems = (fSecondaryPDG == 0 ? fPublication->GetItems() :
    fPublication->GetItemsForSecondary(fSecondaryPDG) );
  TString title = TString::Format("%s(p, %s) reaction at E_{p} = %d MeV",
     gParticlesDAL->GetParticleName(fTargetPDG).Data(),
     gParticlesDAL->GetParticleName(fSecondaryPDG).Data(), en);
  
  gPlotHelper->PrepareCanvas();
  gPlotHelper->ClearCanvas();
  Int_t padsPerRow = gPlotHelper->DivideForNumber(pubItems.size());
  gPlotHelper->DrawBigTitle(title);
  gPlotHelper->SetTitlePosition(33, 0.96, 0.96);
  
  // Load the histograms
  for (UInt_t i = 0; i < fSimulations.size(); ++i)
  {
    G4TData* simulation = gTestDB->LoadData(fSimulations[i], fSecondaryPDG);
    if(simulation == 0 || !simulation->IsLoaded())
    {
      // not loaded, just continue the loop
      continue;
    }
    // prepare histograms and get items
    simulation->BookHistograms(fHnbin, fHlxmin, fHlxmax, 0, i);
    vector<G4TDataItem*> items = (fSecondaryPDG == 0 ? simulation->GetItems() :
      simulation->GetItemsForSecondary(fSecondaryPDG) );
    for (UInt_t j = 0; j < items.size(); ++j)
    {
      G4TDataItem* item = items[j];
      // Set Angle
      Double_t arad = item->GetAngleInRadians(); 
      Double_t dOmega = 2 * fPi * (cos(arad - fDanrad) - cos(arad + fDanrad));
      Double_t w  = (Int_t)simulation->GetCrossSection() * 1000 / dOmega
                     / simulation->GetNumberOfEvents();
      TString Expression = TString::Format("x >> %s", item->GetHistogram()->GetName());
      item->GetData()->Draw(Expression.Data(), "", "goff");
      TH1F* hK = item->GetHistogram();
      hK->Divide(fHDT);
      hK->Add(fHZL);
      hK->Scale(w);
    }
  }
  for(UInt_t i = 0; i < pubItems.size(); ++i)
  {
    G4TDataItem* item = pubItems[i];
    Double_t a  = item->GetCutValue();
    Int_t   mk  = 20;
    // auto size feature
    cout<<"G4TAnalysisTool::InternalExecute: i="<<i<<", padsPerRow = "<<padsPerRow<<endl;
    std::pair<Double_t, Double_t> limits = fPublication->GetLimits(i,padsPerRow);
    if(limits.first && limits.second)
    {
      Double_t rat = std::exp( 0.1 * std::log(limits.first / limits.second) );
      fHymin = limits.first  / rat;
      fHymax = limits.second * rat;
    }

    gPlotHelper->PreparePad(i+1);
    TH1F* frame = gPlotHelper->PrepareFrame(0, fHymin, fHxmax, fHymax,
      TString::Format("%s %g#circ",
        gParticlesDAL->GetParticleName( item->GetSecondaryParticlePDG(), true).Data(), a));
    gPlotHelper->DrawRightAxis(frame, fHxmax, fHymin, fHymax);
    if(!i) // @@ Make a switch
    {
      frame->SetYTitle("d#sigma/pdEd#Omega(mb MeV^{-2} sr^{-1})");
      gPlotHelper->DrawModelsLegend(&fSimulations);
    }
    else if(i == pubItems.size() - 1) frame->SetXTitle("T (MeV)");
    TTree* currentTree  = item->GetData();
    currentTree->Draw(item->GetInvariantFunction().Data(),"","goff");
    TGraphErrors *gr = new TGraphErrors(currentTree->GetSelectedRows(),
                                        currentTree->GetV2(), currentTree->GetV1(),
                                        0,                    currentTree->GetV3() );
    gr->SetMarkerStyle(mk);
    gr->Draw("p");
    for (UInt_t x = 0; x < fSimulations.size(); ++x)
    {
      if(!fSimulations[x]->IsLoaded()) continue; // not loaded, just continue the loop
      G4TDataItem* sim = fSimulations[x]->GetItem(item->fHeader.fSecondaryParticlePDG, a);
      if(!sim)
      {
        TH1F* hist = sim->GetHistogram();
        RenderHSolid(hist, fHfbin, fHnbin,
                     gParticlesDAL->GetParticleMass(item->fHeader.fSecondaryParticlePDG),
                     fSimulations[x]->GetRenderColor(), true);
      }
      else cout<<"Error>G4TAnalysisTool::InternalExecute: Simulation isn't found!"<<endl;
    }
  }
}

//______________________________________________________________________________
void G4TAnalysisTool::SetSecondaryToAnalyze(Int_t idx, Int_t np)
{
  fSecondaryPDG = 0;
  if(idx > 1000000000 || idx < -1000000000) // @@ take into account negative strangeness
  {
    // is PDG code
    vector<Int_t> fragments = fPublication->GetSecondaryPDGs();
    for(UInt_t i = 0; i < fragments.size(); ++i)
    {
      if(fragments[i] == idx)
      {
        fSecondaryPDG = idx;
        //fSecondaryIndex = i - 1;
        break;
      }
    }
    if(!fSecondaryPDG)
    {
      cout<<"Error>G4TAnalysisTool::SetSecondaryToAnalyze: Projectile PDG= "<<idx
          <<" isn't defined!"<<endl;
      return;
    }
  }
  else if(idx > 0)
  {
    // is index
    vector<Int_t> fragments = fPublication->GetSecondaryPDGs();
    if((UInt_t)idx < fragments.size()) fSecondaryPDG = fragments[idx - 1];
    else
    {
      cout<<"Error>G4TAnalysisTool::SetSecondaryToAnalyze: Projectile index# "<<idx
          <<" isn't defined!" << endl;
      return;
    }
  }
  else if(!idx) fSecondaryPDG = 0; // All particles
  else
    cout<<"Error>G4TAnalysisTool::SetSecondaryToAnalyze:UnknownProjectile idx="<<idx<<endl;
  return;
}

