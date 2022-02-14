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
/// \file medical/dna/range/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 72238 2013-07-12 08:40:30Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager() :
    fFileName("AuNP"),
    fpDetector(0)
{
  fpDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());

  Book();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);

  // Define histograms start values
  const G4String nameh1[] = {"h1Events",
                             "h1Edep",
                             "h1SecEnergyNP_charged",
                             "h1SecEnergyNP_nutral",
                             "h1SecEnergyNPSurf_charged",
                             "h1SecEnergyNPSurf_nutral",
                             "h1Sec_charged",
                             "h1Sec_nutral",
                             "h1Chem_0",
                             "h1Chem_1",
                             "h1Chem_2",
                             "h1Chem_3",
                             "h1Chem_4",
                             "h1Chem_5",
                             "h1Chem_6",
                             "h1Chem_7",
                             "h1IncEnergyNPSurf_Back",
                             "h1IncEnergyNPSurf_Fowared"
                            };
  const G4String nameh2[] = {"h2Edep"  ,
                             "h2SecEnergyAbs_charged",
                             "h2SecEnergyAbs_nutral"
                             };

  const G4String titleh1[] = {"Events", 
                              "Energy Deposit Distribution", 
                              "Secondary Energy in NP (charged)",
                              "Secondary Energy in NP (nutral)",
                              "Secondary Energy at NP Surface (charged)",
                              "Secondary Energy at NP Surface (nutral)",
                              "Number of Secondaries in Absorber (charged)",
                              "Number of Secondaries in Absorber (nutral)",
                              "Number of Chemical for ID=0 at 1 psec",
                              "Number of Chemical for ID=1 at 1 psec",
                              "Number of Chemical for ID=2 at 1 psec",
                              "Number of Chemical for ID=3 at 1 psec",
                              "Number of Chemical for ID=4 at 1 psec",
                              "Number of Chemical for ID=5 at 1 psec",
                              "Number of Chemical for ID=6 at 1 psec",
                              "Number of Chemical for ID=7 at 1 psec",
                              "Energy of Incident particle at backwared of GNP",
                              "Energy of Incident particle at forwared of GNP"
                             };
  const G4String titleh2[] = {"Energy Deposit Distribution", 
                              "Secondary Energy vs distance (charged)",
                              "Secondary Energy vs distance (nutral)"
                             };

  // for event counting
  G4int    nbin_eve = 1;
  G4double vmin_eve = 0.;
  G4double vmax_eve = 1.;

  // for SecENP
  G4int    nbin_senp = 1000;
  G4double vmin_senp = 1.;
  G4double vmax_senp = 1000000;

  G4int    NAzm    = fpDetector->GetNReplicaAzm();
  G4int    NR      = fpDetector->GetNReplicaR();
  G4double Rmin    = fpDetector->GetNPRadius() /CLHEP::nm;
  G4double Rmax    = fpDetector->GetAbsRadius()/CLHEP::nm;

  G4int    Runit   = (G4int) (Rmax-Rmin)/NR ;
  NR = NR+(G4int)(Rmin/Runit);

  //for dose distribution
  G4int    nbinAzm    = NAzm;
  G4double vminAzm    = 0. ;
  G4double vmaxAzm    = 360;
  G4int    nbinR2D    = NR;
  G4double vminR2D    = 0.  ;
  G4double vmaxR2D    = 1000.;
  G4int    nbinR      = NR;
  G4double vminR_log  = 10 ;
  G4double vmaxR_log  = Rmax;

  analysisManager->CreateH1(nameh1[ 0], titleh1[ 0], nbin_eve ,
                            vmin_eve , vmax_eve );                   
  analysisManager->CreateH1(nameh1[ 1], titleh1[ 1], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[ 2], titleh1[ 2], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[ 3], titleh1[ 3], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[ 4], titleh1[ 4], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[ 5], titleh1[ 5], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[ 6], titleh1[ 6], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[ 7], titleh1[ 7], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[16], titleh1[16], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[17], titleh1[17], nbin_senp,
                            vmin_senp, vmax_senp,"none","none","log");
  analysisManager->CreateH1(nameh1[ 8], titleh1[ 8], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[ 9], titleh1[ 9], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[10], titleh1[10], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[11], titleh1[11], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[12], titleh1[12], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[13], titleh1[13], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[14], titleh1[14], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");
  analysisManager->CreateH1(nameh1[15], titleh1[15], nbinR    ,
                            vminR_log, vmaxR_log,"none","none","log");

  analysisManager->CreateH2(nameh2[0], titleh2[0], nbinAzm , vminAzm ,
                            vmaxAzm,nbinR2D  ,vminR2D,vmaxR2D);       
  analysisManager->CreateH2(nameh2[1], titleh2[1], nbinR2D , vminR2D ,
                            vmaxR2D,nbin_senp, vmin_senp, vmax_senp,
                            "none","none","none","none","linear","log");
  analysisManager->CreateH2(nameh2[2], titleh2[2], nbinR2D , vminR2D ,
                            vmaxR2D,nbin_senp, vmin_senp, vmax_senp,
                            "none","none","none","none","linear","log");

  analysisManager->SetH1Activation( 0, true);
  analysisManager->SetH1Activation( 1, true);
  analysisManager->SetH1Activation( 2, true);
  analysisManager->SetH1Activation( 3, true);
  analysisManager->SetH1Activation( 4, true);
  analysisManager->SetH1Activation( 5, true);
  analysisManager->SetH1Activation( 6, true);
  analysisManager->SetH1Activation( 7, true);
  analysisManager->SetH1Activation( 8, true);
  analysisManager->SetH1Activation( 9, true);
  analysisManager->SetH1Activation(10, true);
  analysisManager->SetH1Activation(11, true);
  analysisManager->SetH1Activation(12, true);
  analysisManager->SetH1Activation(13, true);
  analysisManager->SetH1Activation(14, true);
  analysisManager->SetH1Activation(15, true);
  analysisManager->SetH1Activation(16, true);
  analysisManager->SetH1Activation(17, true);
  analysisManager->SetH2Activation( 0, true);
  analysisManager->SetH2Activation( 1, true);
  analysisManager->SetH2Activation( 2, true);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
