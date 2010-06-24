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

#include "G4TSimulationTool.h"

G4TSimulationTool *gSimulationTool = new G4TSimulationTool();

ClassImp(G4TSimulationTool)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
Int_t G4TSimulationTool::Run(TString const& publicationFile, Double_t crossSection,
                             const TString& modelName,       const TString& postfix,
                             Int_t runsNumber,               Bool_t useExistingData,
                             TString const& dbPath,          TString const& printQue)
{
  G4TSimHelper::LoadLibraries();

  gBenchmark->Start("Overall Benchmark");
  gROOT->Reset();

  // parse the header and load the publication
  cout<< "G4TSimulationTool::Run: Selected Publication = " << publicationFile <<endl;

  if(fPublication != 0) delete fPublication;
  fPublication = gTestDB->LoadData(publicationFile);

  // prepare arguments
  fProjectilePDG  = fPublication->fHeader.fProjectilePDG;
  fTargetPDG   = fPublication->fHeader.fTargetPDG;
  Double_t pMeV = fPublication->GetProjMomentum();
  cout<<"G4TSimulationTool::Run:"<<fProjectilePDG<<"("<<pMeV<<" MeV/c)+"<<fTargetPDG<<endl;
  fModelName  = modelName;
  fPostfix  = postfix;
  fRunsNumber  = runsNumber;
  fEventsNumber = 2000 * runsNumber;
  fCrossSection = crossSection;
  fUseExistingData= useExistingData;
  gTestDB->SetDirectory(dbPath);

  // get nn and np - Not necessary in the new design
  // Int_t nn = gParticlesDAL->GetN(fTargetPDG);
  // Int_t np = gParticlesDAL->GetZ(fTargetPDG);

  InternalExecute(pMeV, fPublication->fHeader.fSigValue, printQue,2,3, modelName);

  gBenchmark->Show("Overall Benchmark");
  return 0;
}

//______________________________________________________________________________
void G4TSimulationTool::InternalExecute(Double_t mom, Double_t sig, const TString& pq,
                                        Int_t nzone, Int_t nvex, const TString& model,
                                        const TString& dir)
{
  // Prepare  the output simulation @@ Why fModelName cann't be used instead of "model"
  Int_t p = static_cast<Int_t>(mom);
  if(!fSimulation) delete fSimulation;              // *publ*
  fSimulation = new G4TData(fProjectilePDG, fTargetPDG, false, fModelName, fPostfix, p_mom,
                            p, MeV, sig, mBarn, 6);

  fSimulation->SetDirectory(gTestDB->GetDirectory());
  fSimulation->SetCrossSection(fCrossSection);
  fSimulation->SetNumberOfEvents(fEventsNumber);

  Initialize();

  if(!fPublication || !fPublication->IsLoaded())
      cout<<"Warning>G4TSimulationTool::InternalExecute: Publication isn't defined!"<<endl;
  else
  {
    Int_t nofRuns = fRunsNumber;
    // Execute CHIPS
    TString file = gParticlesDAL->GetFileName(fTargetPDG);
    if(fUseExistingData) nofRuns = 0;
    TFile* inFile=fSimHelper.ExecuteTest(fProjectilePDG, fTargetPDG, mom, nofRuns,
                                         fEventsNumber, dir, model);
    Plot( (TTree*)inFile->Get("TuplInclTree"), mom,  fCrossSection, file);
  }
}

//______________________________________________________________________________
void G4TSimulationTool::Plot(TTree* inclTree, Double_t mom, Double_t cs,
                             const TString& file)
{
  if(!fPublication || !fPublication->IsLoaded())
  {
    cout<<"G4TSimulationTool::Plot: Publication is not defined! Aborting..."<<endl;
    return;
  }
  // Get number of events
  Int_t  nevts = fSimHelper.GetEventsNumber();
  // Prepare the histograms
  fHxmax = fPublication->GetMaxT() * 1.03; // Step a bit to the right from the max point
  fHlxmax = std::log(fHxmax);
  cout<<"G4TSimulationTool::Plot: Cuts & histograms, xi="<<fHlxmin<<", xa="<<fHlxmax<<endl;
  PrepareHistograms(fHnbin, fHlxmin, fHlxmax); // Taken out of the G4TTool Constructor
  fPublication->BookHistograms(fHnbin, fHlxmin, fHlxmax, 0);
  cout<<"G4TSimulationTool::Plot: Histograms are prepared"<<endl;
  gPlotHelper->PrepareCanvas();
  cout<<"G4TSimulationTool::Plot: Canvas are prepared"<<endl;
  vector<G4TDataItem*> items = fPublication->GetItems();
  cout<<"G4TSimulationTool::Plot: #ofItems in Publication = "<<items.size()<<endl;
  for(UInt_t x = 0; x < items.size(); ++x)
  {
    G4TDataItem* item = items[x];
    Int_t secondaryPDG = item->GetSecondaryParticlePDG();
    Double_t arad = item->GetAngleInRadians();
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Angle in rad = "<<arad<<endl;
    // Prepare Cut
    TCut cutAtan = TString::Format("abs(atan2(sqrt(Px*Px+Py*Py),Pz)-%g)<%g",
                                   arad, fDanrad).Data();
    TCut cutForSec(gParticlesDAL->GetCut(secondaryPDG));
    TCut cutK = cutForSec && cutAtan;
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Cuts are made = "
        <<gParticlesDAL->GetCut(secondaryPDG)<<endl;
    TString Expression = TString::Format("(log(E-m+0.001)) >> %s",
                                         item->GetHistogram()->GetName() );
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Expression = "<<Expression<<endl;
    inclTree->Draw(Expression.Data(), cutK, "goff");
    //Save Data
    G4TDataItem* simulationItem = fSimulation->AddItem(item->GetHeader());
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Data are saved"<<endl;
    // Copy the tree to the simulation item
    gROOT->cd(); // memory
    TTree* simData = new TTree(item->GetHeader().Data(), item->GetHeader().Data());
    while(!simData) gSystem->Sleep(10);
    Double_t* ptrV1 = inclTree->GetV1();
    Double_t  valV1 = (*ptrV1);
    simData->Branch("x", &valV1,"x/D");
    simData->Fill();
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": SimOUT ASCII are saved to TTree"<<endl;
    for(UInt_t ix = 0; ix < inclTree->GetSelectedRows(); ++ix)
    {
      ptrV1++;
      valV1 = (*ptrV1);
      simData->Fill();
    }
    simulationItem->SetData(simData);
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Simulation Item is set"<<endl;
    // scale @@ It should be modified for different cut values
    Double_t dOmega = 2 * fPi * (cos(arad - fDanrad) - cos(arad + fDanrad));
    Double_t w  = cs / dOmega / nevts;
    cout<<"***>>G4TSimulationTool::Plot: It#"<<x<<": dOm="<<dOmega<<", weight="<<w<<", cs="
        <<cs<<", Nev="<<nevts<<", a="<<arad<<", da="<<fDanrad<<", pi="<<fPi<<endl;
    TH1F* hK = item->GetHistogram();
    // debug print
    Double_t* vW0 = new Double_t[fHnbin];
    Double_t* vW1 = new Double_t[fHnbin];
    Double_t* vW2 = new Double_t[fHnbin];
    Double_t* vW3 = new Double_t[fHnbin];
    Double_t* vW4 = new Double_t[fHnbin];
    for(Int_t i = 0; i < fHnbin; ++i)
    {
      vW0[i] = hK->GetBinCenter(i+1);      // ln(T)
      vW1[i] = hK->GetBinContent(i+1);     // N
    }
    //
    hK->Divide(fHDT);
    // debug print
    for(Int_t i = 0; i < fHnbin; ++i) vW2[i] = hK->GetBinContent(i+1);     // N/dT
    //
    hK->Add(fHZL);
    // debug print
    for(Int_t i = 0; i < fHnbin; ++i) vW3[i] = hK->GetBinContent(i+1);     // N/dT+eps
    //
    hK->Scale(w);
    // debug print
    for(Int_t i = 0; i < fHnbin; ++i) vW4[i] = hK->GetBinContent(i+1);     // Sig/dT/dO
    cout<<"G4TSimulationTool::Plot: Item#"<<x<<": Histogram is modified"<<endl;
    hK->Print(); // temporary print of the histogram
    for(Int_t i = 0; i < fHnbin; ++i) cout<<i<<" : lT="<<vW0[i]<<",N="<<vW1[i]<<", NdT="
                                    <<vW2[i]<<", NdTe="<<vW3[i]<<", dsdTdO="<<vW4[i]<<endl;
    //
    cout<<"G4TSimulationTool::Plot: Yi="<<hK->GetMinimum()<<",Ya="<<hK->GetMaximum()<<endl;
  }
  vector<Int_t> fragments = fPublication->GetSecondaryPDGs();
  cout<<"G4TSimulationTool::Plot: Get Fragments, #ofFrag = "<<fragments.size()<<endl;
  gPlotHelper->ClearCanvas();
  cout<<"G4TSimulationTool::Plot: "<<gParticlesDAL->GetParticleName(fTargetPDG,true).Data()
      <<", proj="<<gParticlesDAL->GetParticleName(fProjectilePDG, true).Data()<<", P="<<mom
      <<", theta="<<gPlotHelper->PrepareCutEnumeration(fPublication).Data()<<", model="
      <<fModelName.Data()<<", Comment="<<fPostfix.Data()<<endl;
  TString tit=TString::Format(" %s(%s,f)X, P=%g MeV/c, #theta=%s#circ (%s_%s)",
                              gParticlesDAL->GetParticleName(fTargetPDG, true).Data(),
                              gParticlesDAL->GetParticleName(fProjectilePDG,true).Data(),
                              mom,
                              gPlotHelper->PrepareCutEnumeration(fPublication).Data(),
                              fModelName.Data(),
                              fPostfix.Data() );
  cout<<"G4TSimulationTool::Plot: Title = "<<tit<<endl;
  gPlotHelper->DrawBigTitle(tit);
  Int_t padsPerRow = gPlotHelper->DivideForNumber(fragments.size());
  cout<<"G4TSimulationTool::Plot: Title is drown, padsPerRow = "<<padsPerRow<<endl;
  gPlotHelper->SetTitlePosition(33, 0.96, 0.96);
  // LOOP over secondary fragments
  for(UInt_t j = 0; j< fragments.size(); ++j)
  {
    gPlotHelper->PreparePad(j+1);
    cout<<"G4TSimulationTool::Plot: Frag#"<<j<<": Pad is prepared"<<endl;
    // auto size feature
    std::pair<Double_t, Double_t> limits = fPublication->GetLimits(j, padsPerRow);
    Double_t cMin=limits.first;
    Double_t cMax=limits.second;
    cout<<"G4TSimulationTool::Plot: Frag#"<<j<<": cMin="<<cMin<<", cMax="<<cMax<<endl;
    if(cMin && cMax)
    {
      // avoid negative or 0 log scale
      if(cMin <= 0) cMin = 0.00000001;
      if(cMax <= 0) cMax = 0.00000002;
      Double_t ria=std::exp( 0.10 * std::log(cMax / cMin) );
      cout<<"G4TSimulationTool::Plot: ria="<<ria<<", cMin="<<cMin<<", cMax="<<cMax<<endl;
      fHymin = cMin / ria;
      fHymax = cMax * ria;
    }
    else cout<<"Warning>G4TSimulationTool::Plot: cMin = cMax = 0"<<endl;
    cout<<"G4TSimulationTool::Plot:ymin="<<fHymin<<",ymax="<<fHymax<<",xma="<<fHxmax<<endl;
    Int_t pdgCode = fragments[j];
    Double_t mass  = gParticlesDAL->GetParticleMass(pdgCode);
    TString Title = gParticlesDAL->GetParticleName(pdgCode);
    cout<<"G4TSimulationTool::Plot:F#"<<j<<":C="<<pdgCode<<",M="<<mass<<",N="<<Title<<endl;
    TH1F* frame = gPlotHelper->PrepareFrame(0,fHymin,fHxmax,fHymax,Title);
    if(j == 0) gPlotHelper->DrawAnglesLegend(fPublication);// Draw legends on the 1st pad
    gPlotHelper->DrawRightAxis(frame, fHxmax, fHymin, fHymax);
    if(!j) frame->SetYTitle("d#sigma/pdEd#Omega(mb MeV^{-2} sr^{-1})");
    if(j == fragments.size() - 1) frame->SetXTitle("T (MeV)");
    gPad->Update();
    // Loop over angles
    vector<G4TDataItem*> sitems = fPublication->GetItemsForSecondary(pdgCode);
    cout<<"G4TSimulationTool::Plot: Axises & Legends are drown, sit="<<sitems.size()<<endl;
    for(UInt_t i = 0; i < sitems.size(); ++i)
    {
      gPad->Update();
      G4TDataItem* item = sitems[i];
      Int_t mk = gPlotHelper->GetMarker(i);
      cout<<"G4TSimulationTool::Plot: Frag#"<<j<<",i#"<<i<<",mar="<<mk<<", m="<<mass<<endl;
      TTree* currentTree  = item->GetData();
      currentTree->Draw(item->GetInvariantFunction().Data(),"","goff");
      TGraphErrors* gr = new TGraphErrors(currentTree->GetSelectedRows(),
                                          currentTree->GetV2(),
                                          currentTree->GetV1(),
                                          0,
                                          currentTree->GetV3());
      gr->SetMarkerStyle(mk);
      gr->Draw("p");
      RenderHSolid( item->GetHistogram(), fHfbin, fHnbin, mass );
    }
    cout<<"G4TSimulationTool::Plot: Frag#"<<j<<": End of the frag loop"<<endl;
  }
  // Save to the DB
  gTestDB->SaveData(fSimulation);
  cout<<"G4TSimulationTool::Plot: End - data are saved"<<endl;
}
