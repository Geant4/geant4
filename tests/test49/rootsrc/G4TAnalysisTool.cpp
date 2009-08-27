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
int G4TAnalysisTool::Run(TString const& publicationFile, Int_t secondaryPDGorIdx,
                         TString const& printQue )
{
 cout << "Loading Libraries..." << endl;
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
void G4TAnalysisTool::InternalExecute(Int_t secondaryPDGorIdx, Int_t np, Int_t nn, Int_t e,
                           const TString& pq, Int_t nzone , Int_t nvex, const TString& dir)
{
  if(fPublication == 0 || !fPublication->IsLoaded())
  {
   cout << "Publication is not defined! Aborting..." << endl;
   return;
  }
  
   SetSecondaryToAnalyze(secondaryPDGorIdx, np);
  
  // Prepare Models
  ArgEnum arg  = fPublication->fHeader.fTypeVar;
  UnitsEnum units = fPublication->fHeader.fTypeUnits;
  Double_t value  = fPublication->fHeader.fTypeValue;
  
  fSimulations.clear();
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "lhep",   arg,
                         value, units,  1));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "chips",  arg,
                         value, units,  2));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "preco",  arg,
                         value, units,  6));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "binary", arg,
                         value, units,  8));
  fSimulations.push_back(new G4TData(fProjectilePDG, fTargetPDG, false, "bertini",arg,
                         value, units,  9));
  
  vector<G4TDataItem*> pubItems = (fSecondaryPDG == 0 ? fPublication->GetItems() :
    fPublication->GetItemsForSecondary(fSecondaryPDG) );
  TString title = TString::Format("%s(p, %s) reaction at E_{p} = %d MeV",
     gParticlesDAL->GetParticleName(fTargetPDG).Data(),
     gParticlesDAL->GetParticleName(fSecondaryPDG).Data(), e);
  
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
    simulation->PrepareHistograms(fHnbin, fHlxmin, fHlxmax, 0, i);
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
  cout << "i = " << i << " pads per row = " << padsPerRow << endl;
  DataItemLimit_t limits = fPublication->GetLimits(i,padsPerRow);
  if(limits.fMin != 0 && limits.fMax != 0)
  {
   fHymin = limits.fMin / exp(  0.10 * log(limits.fMax / limits.fMin) );
   fHymax = limits.fMax * exp(  0.10 * log(limits.fMax / limits.fMin) );
  }

  gPlotHelper->PreparePad(i+1);
  TH1F* frame = gPlotHelper->PrepareFrame(0, fHymin, fHxmax, fHymax,
    TString::Format("%s %g#circ", gParticlesDAL->GetParticleName( item->GetSecondaryParticlePDG(), true).Data(), a));
  gPlotHelper->DrawRightAxis(frame, fHxmax, fHymin, fHymax);

  if(i == 0)
  {
   frame->SetYTitle("d#sigma/pdEd#Omega(mb MeV^{-2} sr^{-1})");
   gPlotHelper->DrawModelsLegend(&fSimulations);
  }
  if(i == pubItems.size() - 1) frame->SetXTitle("T (MeV)");

  TTree* currentTree  = item->GetData();
  currentTree->Draw(item->GetFormula().Data(),"","goff");

  TGraphErrors *gr = new TGraphErrors(currentTree->GetSelectedRows(), currentTree->GetV2(),
    currentTree->GetV1(), 0, currentTree->GetV3());
  gr->SetMarkerStyle(mk);
  gr->Draw("p");

  for (UInt_t x = 0; x < fSimulations.size(); ++x)
  {
   if(!fSimulations[x]->IsLoaded()){
    // not loaded, just continue the loop
    continue;
   }

   G4TDataItem* sim = fSimulations[x]->GetItem(item->fHeader.fSecondaryParticlePDG, a);
   if(sim != 0)
   {
    TH1F* hist = sim->GetHistogram();
    RenderHSolid(hist, fHfbin, fHnbin, gParticlesDAL->GetParticleMass(item->fHeader.fSecondaryParticlePDG), fSimulations[x]->GetRenderColor(),true);
   }else{
    cout << "Error: simulation item not found!" << endl;
   }
  }
 }

      
}


//______________________________________________________________________________
void G4TAnalysisTool::SetSecondaryToAnalyze(Int_t idx, Int_t np)
{
 fSecondaryPDG = 0;
 if(idx > 1000000000 || idx < -1000000000){
  // is PDG code
  vector<Int_t> fragments = fPublication->GetSecondaryPDGs();
  for(UInt_t i = 0; i < fragments.size(); ++i)
  {
   if(fragments[i] == idx){
    fSecondaryPDG = idx;
    //fSecondaryIndex = i - 1;
    break;
   }
  }

  if(fSecondaryPDG == 0)
  {
   cout << "*** Projectile #" << idx << " is not defined!" << endl;
   return;
  }
 }else if(idx > 0){
  // is index
  vector<Int_t> fragments = fPublication->GetSecondaryPDGs();
  if((UInt_t)idx < fragments.size()){
   fSecondaryPDG = fragments[idx - 1];
  }else{
   cout << "*** Projectile #" << idx << " is not defined!" << endl;
   return;
  }
 }
 else if(idx == 0)
 {
  // All particles
  fSecondaryPDG = 0;
 }
 else
 {
  cout << "*** Projectile #" << idx << " is not defined!" << endl;
 }

}

