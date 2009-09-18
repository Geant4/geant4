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
// $Id: RunAction.cc,v 1.8 2009-09-18 17:34:54 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4Geantino.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim,
                     HistoManager* hist)
:detector(det), primary(prim), histoManager(hist)
{  
  writeFile = false; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  //
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  //initialize cumulative quantities
  //
  G4int nbPixels = detector->GetSizeVectorPixels();
  G4int size = totalEnergy.size();
  if (size < nbPixels) {
    visibleEnergy.resize(nbPixels);  
      totalEnergy.resize(nbPixels);     

    visibleEnergy2.resize(nbPixels);
      totalEnergy2.resize(nbPixels);
  }
  
  for (G4int k=0; k<nbPixels; k++) {
   visibleEnergy[k] = visibleEnergy2[k] = totalEnergy[k]= totalEnergy2[k] = 0.0;
  }
        
  G4int n1pxl = detector->GetN1Pixels();
  size = layerEtot.size();
  if (size < n1pxl) {  
    layerEvis.resize(n1pxl);
    layerEtot.resize(n1pxl);
    layerEvis2.resize(n1pxl);
    layerEtot2.resize(n1pxl);
  }
  
  for (G4int k=0; k<n1pxl; k++) {
   layerEvis[k] = layerEvis2[k] = layerEtot[k]= layerEtot2[k] = 0.0;
  }
     
  nbEvents = 0;  
  calorEvis = calorEvis2 = calorEtot = calorEtot2 = Eleak = Eleak2 = 0.;
  EdLeak[0] = EdLeak[1] = EdLeak[2] = 0.;
  nbRadLen  = nbRadLen2 = 0.;
                    
  //histograms
  //
  histoManager->book();
  
  //create ascii file for pixels
  //
  if (writeFile) CreateFilePixels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent_1(G4int pixel, G4double Evis, G4double Etot)
{
  //accumulate statistic per pixel
  //
  visibleEnergy[pixel] += Evis;  visibleEnergy2[pixel] += Evis*Evis;
    totalEnergy[pixel] += Etot;    totalEnergy2[pixel] += Etot*Etot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent_2(G4int layer, G4double Evis, G4double Etot)
{
  //accumulate statistic per layer
  //
  layerEvis[layer] += Evis;  layerEvis2[layer] += Evis*Evis;
  layerEtot[layer] += Etot;  layerEtot2[layer] += Etot*Etot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent_3(G4double calEvis, G4double calEtot,
                               G4double eleak)
{
  //accumulate statistic
  //
  nbEvents++;
  calorEvis += calEvis;  calorEvis2 += calEvis*calEvis;
  calorEtot += calEtot;  calorEtot2 += calEtot*calEtot;  
  Eleak += eleak;  Eleak2 += eleak*eleak;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillDetailedLeakage(G4int icase, G4double energy)
{
  //forward, backward, lateral leakage
  //
  EdLeak[icase] += energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillNbRadLen(G4double dn)
{
  //total number of radiation length
  //
  nbRadLen += dn; nbRadLen2 += dn*dn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run*)
{ 
  //calorimeter
  //
  detector->PrintCalorParameters();
      
  //run conditions
  //   
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();

  G4int prec = G4cout.precision(3);
  
  G4cout << " The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") << " through the calorimeter" << G4endl;
	 
  G4cout << "------------------------------------------------------------"
         << G4endl;

  //if no events, return
  //
  if (nbEvents == 0) return;
  
  //compute and print statistic
  //
  std::ios::fmtflags mode = G4cout.flags();
  
  //number of radiation length
  //
  if (particle == G4Geantino::Geantino() ) {
    G4double meanNbRadL  = nbRadLen/ nbEvents;
    G4double meanNbRadL2 = nbRadLen2/nbEvents;
    G4double varNbRadL = meanNbRadL2 - meanNbRadL*meanNbRadL;
    G4double rmsNbRadL = 0.;
    if (varNbRadL > 0.) rmsNbRadL = std::sqrt(varNbRadL);
    G4double effRadL = (detector->GetCalorThickness())/meanNbRadL;
    G4cout.precision(5);
    G4cout
      << "\n Calor : mean number of Rad Length = " 
      << meanNbRadL << " +- "<< rmsNbRadL
      << "  --> Effective Rad Length = "
      << G4BestUnit( effRadL,"Length") << G4endl;    
    
    G4cout << "------------------------------------------------------------"
           << G4endl;        
  }
  
  G4cout.precision(prec);	 
  G4cout << "\n             " 
         << "visible Energy          (rms/mean)        "
         << "total Energy           (rms/mean)" << G4endl;

  G4double meanEvis,meanEvis2,varianceEvis,rmsEvis,resEvis;
  G4double meanEtot,meanEtot2,varianceEtot,rmsEtot,resEtot;
  
  G4int n1pxl = detector->GetN1Pixels();
    
  for (G4int i1=0; i1<n1pxl; i1++) {
    //visible energy
    meanEvis  = layerEvis[i1] /nbEvents;
    meanEvis2 = layerEvis2[i1]/nbEvents;    
    varianceEvis = meanEvis2 - meanEvis*meanEvis;
    rmsEvis = 0.;
    if (varianceEvis > 0.) rmsEvis = std::sqrt(varianceEvis);
    resEvis = 100*rmsEvis/meanEvis;
    histoManager->FillHisto(3, i1+0.5, meanEvis);
         
    //total energy
    meanEtot  = layerEtot[i1] /nbEvents;
    meanEtot2 = layerEtot2[i1]/nbEvents;    
    varianceEtot = meanEtot2 - meanEtot*meanEtot;
    rmsEtot = 0.;
    if (varianceEtot > 0.) rmsEtot = std::sqrt(varianceEtot);
    resEtot = 100*rmsEtot/meanEtot;
    histoManager->FillHisto(4, i1+0.5, meanEtot);    

    //print
    //
    G4cout
      << "\n   layer " << i1 << ": "
      << std::setprecision(5)
      << std::setw(6) << G4BestUnit(meanEvis,"Energy") << " +- "
      << std::setprecision(4)
      << std::setw(5) << G4BestUnit( rmsEvis,"Energy") << "  ("
      << std::setprecision(2) 
      << std::setw(3) << resEvis  << " %)" 
      << "     "
      << std::setprecision(5)
      << std::setw(6) << G4BestUnit(meanEtot,"Energy") << " +- "
      << std::setprecision(4)
      << std::setw(5) << G4BestUnit( rmsEtot,"Energy") << "  ("
      << std::setprecision(2) 
      << std::setw(3) << resEtot  << " %)"; 
  }
  G4cout << G4endl;
  
  //calorimeter: visible energy
  meanEvis  = calorEvis /nbEvents;
  meanEvis2 = calorEvis2/nbEvents;
  varianceEvis = meanEvis2 - meanEvis*meanEvis;
  rmsEvis = 0.;
  if (varianceEvis > 0.) rmsEvis = std::sqrt(varianceEvis);
  resEvis = 100*rmsEvis/meanEvis;
  
  //calorimeter: total energy
  meanEtot  = calorEtot /nbEvents;
  meanEtot2 = calorEtot2/nbEvents;
  varianceEtot = meanEtot2 - meanEtot*meanEtot;
  rmsEtot = 0.;
  if (varianceEtot > 0.) rmsEtot = std::sqrt(varianceEtot);
  resEtot = 100*rmsEtot/meanEtot;
    
  //print
  //
  G4cout
    << "\n   total calor : "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEvis,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEvis,"Energy") << "  ("
    << std::setprecision(2) 
    << std::setw(3) << resEvis  << " %)" 
    << "     "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEtot,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEtot,"Energy") << "  ("
    << std::setprecision(2) 
    << std::setw(3) << resEtot  << " %)";
                     
  G4cout << "\n------------------------------------------------------------"
         << G4endl;

  //leakage
  G4double meanEleak,meanEleak2,varianceEleak,rmsEleak,ratio;
  meanEleak  = Eleak /nbEvents;
  meanEleak2 = Eleak2/nbEvents;
  varianceEleak = meanEleak2 - meanEleak*meanEleak;
  rmsEleak = 0.;
  if (varianceEleak > 0.) rmsEleak = std::sqrt(varianceEleak);
  ratio = 100*meanEleak/energy;
  
  G4double forward = 100*EdLeak[0]/(nbEvents*energy);
  G4double bakward = 100*EdLeak[1]/(nbEvents*energy);
  G4double lateral = 100*EdLeak[2]/(nbEvents*energy);      
  //print
  //
  G4cout
    << "\n   Leakage : "
    << std::setprecision(5)
    << std::setw(6) << G4BestUnit(meanEleak,"Energy") << " +- "
    << std::setprecision(4)
    << std::setw(5) << G4BestUnit( rmsEleak,"Energy") 
    << "\n   Eleak/Ebeam ="
    << std::setprecision(3) 
    << std::setw(4) << ratio  << " %  ( forward ="
    << std::setw(4) << forward  << " %;   backward ="
    << std::setw(4) << bakward  << " %;   lateral ="
    << std::setw(4) << lateral  << " %)"             
    << G4endl;

  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  //save histograms   
  histoManager->save();

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
#include <fstream>

void RunAction::CreateFilePixels()
{
  //create file and write run header
  //
  G4String name = histoManager->GetFileName(); 
  G4String fileName = name + ".pixels.ascii";
  
  std::ofstream File(fileName, std::ios::out);

  G4int n1pxl   = detector->GetN1Pixels();
  G4int n2pxl   = detector->GetN2Pixels();
  G4int n1shift = detector->GetN1Shift();    
  G4int nbEvents    = G4RunManager::GetRunManager()->GetCurrentRun()
                     ->GetNumberOfEventToBeProcessed();
  File << nbEvents << " " << n1pxl << " " <<  n2pxl << " " << n1shift
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


