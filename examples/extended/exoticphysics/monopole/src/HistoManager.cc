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
// $Id: HistoManager.cc 104174 2017-05-15 12:12:45Z selles $
// GEANT4 tag $Name: geant4-09-04-cand-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "DetectorConstruction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager(DetectorConstruction* det, G4double binLength)
  : fDetector(det), fBinLength(binLength)

{
  fVerbose = true;

  fHistoId.resize(MaxHisto);
  fExist.resize(MaxHisto);
  fActive.resize(MaxHisto);
  fLabel.resize(MaxHisto);
  fTitle.resize(MaxHisto);
  fNbins.resize(MaxHisto);
  fVmin.resize(MaxHisto);
  fVmax.resize(MaxHisto);
  fUnit1.resize(MaxHisto);
  fUnit2.resize(MaxHisto);
  fWidth.resize(MaxHisto);
  fIds.resize(MaxHisto);

  // histograms
  fNbHisto = 0;
  for (G4int k=0; k<MaxHisto; k++) {
    fHistoId[k] = 0;
    fExist[k] = false;
    fActive[k] = false;
    fUnit1[k]  = 1.0;
    fUnit2[k]  = 1.0;
    fWidth[k] = 1.0;
  }

  fNtupleActive = false;
  fTupleName  = "tuple";
  fTupleTitle = "test";

  fTupleI.assign(MaxHisto,-1);
  fTupleF.assign(MaxHisto,-1);
  fTupleD.assign(MaxHisto,-1);

  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create or get analysis manager
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(false);  // enable inactivation of histograms

  // Define the histogramms 
  if(analysisManager->GetNofH1s()==0)
    {
      G4double length  = fDetector->GetAbsorSizeX();
      G4int nbBins = G4lrint(length / fBinLength);

      // Create histograms
      fNbHisto = 0;
      Add1D(0,"dummy", nbBins, 0, length, "mm");
      Add1D(1,"Edep (MeV/mm) along absorber (mm)", nbBins, 0, length, "mm");
      Add1D(2,"DEDX (MeV/mm) of proton", 100, -3., 7.);
      Add1D(3,"DEDX (MeV/mm) of monopole", 100, -3., 7.);
      Add1D(4,"Range(mm) of proton", 100, -3., 7., "mm");
      Add1D(5,"Range(mm) of monopole", 100, -3., 7., "mm");

     // Creating an 1-dimensional histograms in the root directory of the tree
      for(G4int i=0; i<fNbHisto; i++) 
          {
            fHistoId[i] = analysisManager->CreateH1(fIds[i], fTitle[i], 
                                                    fNbins[i], fVmin[i], fVmax[i]);
            analysisManager->SetH1Activation(fHistoId[i],fActive[i]);
          }

      // Creating a tuple factory, whose tuples will be handled by the tree
      if(fNtupleActive) {
        analysisManager->CreateNtuple(fTupleName,fTupleTitle); 
        G4int i;
        G4int n = fNtupleI.size();
        for(i=0; i<n; ++i) { 
          if(fTupleI[i] == -1) 
             {  fTupleI[i] = analysisManager->CreateNtupleIColumn(fNtupleI[i]); }
        }
        n = fNtupleF.size();
        for(i=0; i<n; ++i) { 
          if(fTupleF[i] == -1) 
             {  fTupleF[i] = analysisManager->CreateNtupleFColumn(fNtupleF[i]); }
        }
        n = fNtupleD.size();
        for(i=0; i<n; ++i) { 
          if(fTupleD[i] == -1) 
             {  fTupleD[i] = analysisManager->CreateNtupleDColumn(fNtupleD[i]); }
        }
      }
    }

  // Added to catch the SetActivation parameters set througj UI interface
  for (G4int k=0; k<MaxHisto; k++) 
          fExist[k] = analysisManager->GetH1Activation(fHistoId[k]);

   // Check if a filename is set
  if(analysisManager->GetFileName()=="") return;
  
  // Activate the analysisManage ronly if a filename is set
  analysisManager->SetActivation(true);     // enable inactivation of histograms

}

void HistoManager::SetBinLength(G4double binLength)
{
  fBinLength = binLength;

  G4double length  = fDetector->GetAbsorSizeX();
  G4int nbBins = G4lrint(length / fBinLength);
  SetHisto1D(1, nbBins, 0., length, "mm");
}


void HistoManager::Update(G4double binLength)
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4double histo1_binLength = analysisManager->GetH1Width(fHistoId[1]);
  if(fabs(histo1_binLength-binLength)>0.01)
    {
      G4double length  = fDetector->GetAbsorSizeX();
      G4int nbBins = G4lrint(length / binLength);
      analysisManager->SetH1(fHistoId[1], nbBins,
                       G4AnalysisManager::Instance()->GetH1Xmin(fHistoId[1]),
                       G4AnalysisManager::Instance()->GetH1Xmax(fHistoId[1]));
    }                     
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Add1D(G4int histoId, const G4String& name, G4int nb, 
                         G4double x1, G4double x2, 
                         const G4String& u1, const G4String& u2)
{
  std::stringstream sg;
  sg<<histoId;
  const G4String id = sg.str();

  if(fVerbose) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name 
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u1<<" "<<u2 
           << G4endl;
  }
  fNbHisto++;

  double vUnit1= (u1=="none")? 1. : (G4UnitDefinition::GetValueOf(u1));
  fUnit1[histoId]=vUnit1;
  double vUnit2= (u2=="none")? 1. : (G4UnitDefinition::GetValueOf(u2));
  fUnit2[histoId]=vUnit2;

  x1 /= vUnit1;
  x2 /= vUnit1;
  fNbins[histoId]=nb;
  fVmin[histoId]=x1;
  fVmax[histoId]=x2;
  fWidth[histoId] = (x2-x1)/nb;
  fTitle[histoId]=name;
  fIds[histoId]=id;
  G4int exist = (name.substr(0,5)=="dummy"||name.substr(0,5)=="Dummy")?false:true;
  fExist[histoId]=exist;
  fActive[histoId]=exist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto1D(G4int ih,
                 G4int nbins, G4double valmin, G4double valmax, 
                 const G4String& unit)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }

//   const G4String id[] = { "0", "1", "2", "3", "4", "5", "6" };
//   const G4String title[] =
//                 { "dummy",                                        //0
//                   "continuous energy loss along primary track",        //1
//                   "energy from secondaries",                        //2
//                   "total energy lost by primary track",                //3
//                   "energy spectrum of e-+",                        //4
//                   "energy spectrum of gamma",                        //5
//                   "step size"                                        //6
//                  };

//   G4String titl = fTitle[ih];
  G4double vmin = valmin, vmax = valmax;

  if (unit != "none") {
    fUnit1[ih] = G4UnitDefinition::GetValueOf(unit);
  }
  vmin = valmin/fUnit1[ih]; 
  vmax = valmax/fUnit1[ih];

  fExist[ih] = true;
  fNbins[ih] = nbins;
  fVmin[ih]  = vmin;
  fVmax[ih]  = vmax;
  fWidth[ih] = (valmax-valmin)/nbins;

  if(fTitle[ih].substr(0,5)=="dummy") fExist[ih]=false;
  if(fTitle[ih].substr(0,5)=="Dummy") fExist[ih]=false;
  fActive[ih]=fExist[ih];

  G4cout << "----> SetHisto " << ih << ": " << fTitle[ih] << ";  "
         << nbins << " bins from "
         << vmin << " " << unit << " to " << vmax << " " << fUnit1[ih] << " - "
         <<fExist[ih]<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double e, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does not exist; e= " << e << " w= " << weight << G4endl;
    return;
  }

  if(!fExist[ih]) return;
  G4AnalysisManager::Instance()->FillH1( fHistoId[ih], e/fUnit1[ih], weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Scale(G4int ih, G4double fac)
{
 if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::Scale() : histo " << ih
           << "does not exist.  (fac = " << fac << ")" << G4endl;
    return;
  }

  if(!fExist[ih]) return;
  G4AnalysisManager::Instance()->GetH1(fHistoId[ih])->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

