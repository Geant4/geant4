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
//    *******************************
//    *                             *
//    *    RemSimAnalysisManager.cc *
//    *                             *
//    *******************************
//
// $Id: RemSimAnalysisManager.cc,v 1.14 2010-11-10 00:18:04 asaim Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifdef  G4ANALYSIS_USE 
#include <stdlib.h>
#include <fstream>
#include "RemSimAnalysisManager.hh"
#include "G4ios.hh"
#include <AIDA/AIDA.h>
#include "G4RunManager.hh"
#include <vector>

RemSimAnalysisManager* RemSimAnalysisManager::instance = 0;

RemSimAnalysisManager::RemSimAnalysisManager() 
  :  aFact(0), treeFact(0), theTree(0), dataPointFactory(0),
     histogramFactory(0), dataPoint(0), energyDeposit(0),
     primary(0), secondaryDeposit(0),
     primaryInitialE(0), primaryInitialEout(0), initialE(0), 
     initialEout(0), shape(0), energyShape(0) , histo_secondary_phantom(0),histo_secondary(0),
     histo_proton(0), histo_neutron(0), histo_pion(0), histo_alpha(0), 
     histo_positron(0), histo_electron(0), histo_gamma(0), histo_muon(0), histo_other(0), histo_neutrino(),
     histo_proton_reaching(0), histo_neutron_reaching(0), histo_pion_reaching(0), histo_alpha_reaching(0), 
     histo_positron_reaching(0), histo_electron_reaching(0), histo_gamma_reaching(0), histo_muon_reaching(0), 
     histo_other_reaching(0),
     histo_proton_slice(0), histo_neutron_slice(0), histo_pion_slice(0), histo_alpha_slice(0), 
     histo_positron_slice(0), histo_electron_slice(0), histo_gamma_slice(0), 
     histo_muon_slice(0), histo_other_slice(0), histo_vehicle(0)
{ 
  aFact = AIDA_createAnalysisFactory();
  fileFormat = "xml";
}

RemSimAnalysisManager::~RemSimAnalysisManager() 
{ 
  delete histo_vehicle; histo_vehicle = 0;
  delete histo_proton_reaching; histo_proton_reaching = 0;
  delete histo_neutron_reaching; histo_neutron_reaching = 0;
  delete histo_pion_reaching; histo_pion_reaching = 0; 
  delete histo_alpha_reaching; histo_alpha_reaching = 0;
  delete histo_positron_reaching; histo_positron_reaching = 0;
  delete histo_electron_reaching; histo_electron_reaching = 0;
  delete histo_gamma_reaching; histo_gamma_reaching = 0;
  delete histo_muon_reaching; histo_muon_reaching = 0;
  delete histo_other_reaching; histo_other_reaching = 0;
  delete histo_proton_slice; histo_proton_slice = 0;
  delete histo_neutron_slice; histo_neutron_slice = 0;
  delete histo_pion_slice; histo_pion_slice = 0; 
  delete histo_alpha_slice; histo_alpha_slice = 0;
  delete histo_positron_slice; histo_positron_slice = 0;
  delete histo_electron_slice; histo_electron_slice = 0;
  delete histo_gamma_slice; histo_gamma_slice = 0;
  delete histo_muon_slice; histo_muon_slice = 0;
  delete histo_other_slice; histo_other_slice = 0;
  delete histo_proton; histo_proton = 0;
  delete histo_neutron; histo_neutron = 0;
  delete histo_pion; histo_pion = 0; 
  delete histo_alpha; histo_alpha =0;
  delete histo_positron; histo_positron = 0;
  delete histo_electron; histo_electron = 0;
  delete histo_gamma; histo_gamma = 0;
  delete histo_muon; histo_muon = 0;
  delete histo_other; histo_other = 0;
  delete histo_neutrino; histo_neutrino = 0;
  delete histo_secondary; histo_secondary = 0;
  delete histo_secondary_phantom; histo_secondary_phantom = 0;
  delete energyShape; energyShape = 0;
  delete shape; shape = 0;
  delete initialEout; initialEout = 0;
  delete initialE; initialE = 0;
  delete primaryInitialEout; primaryInitialEout = 0; 
  delete primaryInitialE; primaryInitialE = 0;
  delete secondaryDeposit; secondaryDeposit = 0;
  delete primary; primary = 0;
  delete energyDeposit; energyDeposit = 0;
  delete dataPoint; dataPoint = 0; 
  delete histogramFactory; histogramFactory = 0;
  delete dataPointFactory; dataPointFactory = 0;
  delete theTree; theTree = 0;
  delete treeFact; treeFact = 0;
  delete aFact; aFact = 0;
}

RemSimAnalysisManager* RemSimAnalysisManager::getInstance()
{
  if (instance == 0) instance = new RemSimAnalysisManager;
  return instance;
}

void RemSimAnalysisManager::book() 
{ 
  treeFact = aFact -> createTreeFactory();
   

  theTree = treeFact -> create("remsim.xml", "xml",false, true,"");

  G4cout << "The format of the output file is xml" << G4endl;

  histogramFactory = aFact -> createHistogramFactory(*theTree);
  energyDeposit = histogramFactory -> createHistogram1D("10",
                                                        "Energy Deposit",
                                                        30,// number of bins
					                0.,//xmin
                                                        30.);//xmax 
 
  /// A variable binning is created for some histograms
  G4int size  = 1901;

  std::vector<double> vector_bin; 

  G4double x_value = 0.;
  for (G4int ptn = 0; ptn < size; ptn++ ) 
    {
      G4int binNb = 1000;
      if (ptn < binNb)
	{
	  vector_bin.push_back(x_value);
	  x_value += 10; 
	}
      else
	{ vector_bin.push_back(x_value); x_value += 100.;} 
    }

  primary = histogramFactory -> createHistogram1D("20",
						  "Initial energy of primary particles", 
						  vector_bin);
 
  secondaryDeposit = histogramFactory -> createHistogram1D("30",
							   "EnergyDeposit given by secondaries", 
							   30,0.,30.);

  primaryInitialE = histogramFactory -> createHistogram1D("40",
							  "Initial energy of primaries reaching the phantom", vector_bin);

  primaryInitialEout = histogramFactory -> createHistogram1D("50",
							     "Initial energy of primaries ougoing the phantom", 
							     vector_bin);

  initialE = histogramFactory -> createHistogram1D("60",
						   "Energy of primaries reaching the phantom", 
						   vector_bin);

  initialEout = histogramFactory -> createHistogram1D("70",
						      "Energy of primaries outgoing the phantom", 
						      vector_bin);

  shape =  histogramFactory -> createHistogram2D("80",
						 "Shape", 
						 320,-160.,160.,
						 320, -160.,160.);

  energyShape = histogramFactory -> createHistogram2D("90", 
						      "energyDepShape",
						      320, -160.,160.,
						      320, -160.,160.);

  histo_secondary_phantom = histogramFactory -> createHistogram1D("100", 
								  "secondary particles produced in the phantom",
								  10, 0., 10.);

  histo_secondary = histogramFactory -> createHistogram1D("300", 
							  "secondary particles reaching the phantom",
							  10, 0., 10.);

  histo_vehicle = histogramFactory -> createHistogram1D("500", 
                                                        "secondary particles in the vehicle",
                                                        10, 0., 10.);

  histo_proton = histogramFactory -> createHistogram1D("110","Energy of secondary p produced in the phantom", 
						       100, 0., 1000.); // up to 1 GeV, 10 MeV bin

  histo_neutron= histogramFactory -> createHistogram1D("120","Energy of secondary n produced in the phantom", 
						       100, 0., 1000.); // up to 1 GeV, 10 MeV bin

  histo_pion = histogramFactory-> createHistogram1D("130","Energy of secondary pions produced in the phantom",
						    200, 0., 2000.);// up to 2 GeV, 10 MeV bin


  histo_alpha = histogramFactory-> createHistogram1D("140","Energy of secondary alpha produced in the phantom", 
						     100, 0.,100.); // up to 100 MeV, 1 MeV bin

  histo_positron = histogramFactory-> createHistogram1D("150","Energy of secondary e+ produced in the phantom", 
							100, 0., 1000.); // up to 1 GeV, 10 MeV bin

  histo_electron = histogramFactory-> createHistogram1D("160","Energy of secondary e- produced in the phantom", 
							1000, 0., 10.); // up to 10 MeV, 10. keV bin

  histo_gamma = histogramFactory-> createHistogram1D("170","Energy of secondary gamma produced in the phantom",
						     100, 0., 10.); // up to 10 MeV, 1 MeV bin

  histo_muon = histogramFactory-> createHistogram1D("180","Energy of secondary mu produced in the phantom", 
						    100, 0., 1000.); // up to 1 GeV, 10 MeV bin

  histo_other= histogramFactory -> createHistogram1D("190","Energy of secondary other particles produced in the phantom", 
						     200,0., 200.); // up to 200 MeV, 1 MeV bin

  histo_neutrino = histogramFactory -> createHistogram1D("400", "Energy of secondary neutrinos produced in the phantom", 1000, 0., 10000.);

  histo_proton_slice = histogramFactory -> createHistogram1D("200",
							     " Phantom Slice where  secondary protons are produced", 
							     30, 0., 30.);

  histo_neutron_slice = histogramFactory -> createHistogram1D("210",
							      "Phantom Slice where  secondary neutrons are produced", 
							      30, 0., 30.);

  histo_pion_slice = histogramFactory -> createHistogram1D("220",
							   "Phantom Slice where  secondary pions are produced", 
							   30, 0., 30.);

  histo_alpha_slice = histogramFactory -> createHistogram1D("230",
							    "Phantom Slice where  secondary alpha are produced",
							    30, 0., 30.);

  histo_positron_slice = histogramFactory -> createHistogram1D("240",
							       "Phantom Slice where  secondary positrons are produced ",
							       30, 0.,30.);

  histo_electron_slice = histogramFactory -> createHistogram1D("250",
							       "Phantom Slice where  secondary electrons are produced", 
							       30, 0., 30.);

  histo_gamma_slice = histogramFactory -> createHistogram1D("260",
							    "Phantom Slice where  secondary gamma are produced", 
							    30, 0., 30.);

  histo_muon_slice = histogramFactory -> createHistogram1D("270",
							   "Phantom Slice where  secondary muons are produced", 
							   30, 0., 30.);

  histo_other_slice = histogramFactory -> createHistogram1D("280",
							    " Phantom Slice where  secondary other particles are produced",
							    30,0.,30.);

  histo_proton_reaching = histogramFactory -> createHistogram1D("310","Energy of secondary p reaching the phantom", 
								1000, 0., 10000.); // up to 10 GeV, 10 MeV bin

  histo_neutron_reaching= histogramFactory -> createHistogram1D("320","Energy of secondary n reaching the phantom", 
								1000, 0., 10000.); // up to 10 GeV, 10 MeV bin

  histo_pion_reaching = histogramFactory -> createHistogram1D("330","Energy of secondary pions reaching in the phantom", 
							     1000, 0., 10000.); // up to 10 GeV, 10 MeV bin

  histo_alpha_reaching = histogramFactory -> createHistogram1D("340","Energy of secondary alpha reaching the phantom", 
							      100, 0., 100.); // up to 100 MeV, 1 MeV bin

  histo_positron_reaching = histogramFactory -> createHistogram1D("350","Energy of secondary e+ reaching the phantom", 
								 500, 0., 5000.); // up to 5 GeV, 10 MeV bin

  histo_electron_reaching = histogramFactory -> createHistogram1D("360","Energy of secondary e- reaching the phantom",
								 200, 0., 2000.); //up to 2 GeV, 10 MeV bin 

  histo_gamma_reaching = histogramFactory -> createHistogram1D("370","Energy of secondary gamma reaching  the phantom",
							      200, 0., 2000.);//up to 2 GeV, 10 MeV bin

  histo_muon_reaching = histogramFactory -> createHistogram1D("380","Energy of secondary mu reaching the phantom", 
							     200, 0., 2000.); // up to 2 GeV, 10 MeV bin

  histo_other_reaching= histogramFactory -> createHistogram1D("390","Energy of secondary other particles reaching the phantom", 
							      100, 0., 1000.); // up to 1 GeV, 10 MeV bin
}

void RemSimAnalysisManager::energyDepositStore(G4int layer, G4double eDep)
{
  energyDeposit -> fill(layer,eDep);
}

void RemSimAnalysisManager::primaryParticleEnergyDistribution(G4double energy)
{
  primary -> fill(energy);
}

void RemSimAnalysisManager::SecondaryEnergyDeposit(G4int i, G4double EnergyDep)
{
  secondaryDeposit -> fill(i,EnergyDep); 
}

void RemSimAnalysisManager::PrimaryInitialEnergyIn(G4double EnergyDep)
{
  primaryInitialE -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryInitialEnergyOut(G4double EnergyDep)
{
  primaryInitialEout -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryEnergyIn(G4double EnergyDep)
{
  initialE -> fill(EnergyDep); 
}

void RemSimAnalysisManager::PrimaryEnergyOut(G4double EnergyDep)
{
  initialEout -> fill(EnergyDep); 
}

void RemSimAnalysisManager::particleShape(G4double x, G4double y)
{
  shape -> fill(x,y); 
}

void RemSimAnalysisManager::energyDepShape(G4double x, G4double y, G4double energyDep)
{
  energyShape -> fill(x,y, energyDep); 
}

void RemSimAnalysisManager:: SetFormat(G4String format)
{ 
  //  fileFormat = format;

  if (fileFormat != "hbook")
  G4cout << format << "is not available at the moment !!!" << G4endl; 
}
void RemSimAnalysisManager:: SecondaryInPhantom(G4int i)
{
 histo_secondary_phantom -> fill(i);

 // if i = 0  secondary proton, i = 1 neutron, i= 2 pion, i= 3 alpha, 
 // i =4 other, i=5 electron, i = 6 gamma, i= 7 positrons, 
 // i = 8 muons, i= 9 neutrinos
}
 
void  RemSimAnalysisManager::SecondaryProtonInPhantom(G4double energy)
{
 histo_proton -> fill(energy);
 //G4cout<< "proton energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryNeutronInPhantom(G4double energy)
{ 
  histo_neutron -> fill(energy);
  //G4cout<< "neutron energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryPionInPhantom(G4double energy)
{
  histo_pion -> fill(energy);
 //G4cout<< "pion energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryAlphaInPhantom(G4double energy)
{
  histo_alpha -> fill(energy);
  //G4cout<< "alpha energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryPositronInPhantom(G4double energy)
{
  histo_positron -> fill(energy);
  //G4cout<< "positron energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryElectronInPhantom(G4double energy)
{
  histo_electron -> fill(energy);
  //G4cout<< "electron energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryGammaInPhantom(G4double energy)
{ 
  histo_gamma -> fill(energy);
  //G4cout<< "gamma energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryMuonInPhantom(G4double energy)
{ 
  histo_muon -> fill(energy);
  //G4cout<< "muon energy : "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryOtherInPhantom(G4double energy)
{ 
  histo_other -> fill(energy);
  //G4cout<< "other energy : "<< energy << G4endl;
}
void RemSimAnalysisManager::SecondaryNeutrinoInPhantom (G4double energy)
{
  histo_neutrino -> fill(energy);
  // G4cout<< "neutrino energy:" << energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryProtonInPhantomSlice(G4double slice)
{
  histo_proton_slice -> fill(slice);
 //G4cout<< "proton slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryNeutronInPhantomSlice(G4double slice)
{ 
  histo_neutron_slice -> fill(slice);
  //G4cout<< "neutron slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryPionInPhantomSlice(G4double slice)
{
  histo_pion_slice -> fill(slice);
  //G4cout<< "pion slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryAlphaInPhantomSlice(G4double slice)
{
  histo_alpha_slice -> fill(slice);
  //G4cout<< "alpha slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryPositronInPhantomSlice(G4double slice)
{
  histo_positron_slice -> fill(slice);
  //G4cout<< "positron slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryElectronInPhantomSlice(G4double slice)
{
  histo_electron_slice -> fill(slice);
  //G4cout<< "electron slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryGammaInPhantomSlice(G4double slice)
{ 
 histo_gamma_slice -> fill(slice);
 //G4cout<< "gamma slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryMuonInPhantomSlice(G4double slice)
{ 
  histo_muon_slice -> fill(slice);
  //G4cout<< "muon slice : "<< slice << G4endl;
}

void  RemSimAnalysisManager::SecondaryOtherInPhantomSlice(G4double slice)
{ 
  histo_other_slice -> fill(slice);
  //G4cout<< "other slice : "<< slice << G4endl;
}

void RemSimAnalysisManager::SecondaryReachingThePhantom(G4int i)
{
 histo_secondary -> fill(i); 
 // if i =0  secondary proton, i =1 neutron, i=2 pion, i=3 alpha, 
 // i =4 other, i=5 electron, i = 6 gamma, i=7 positrons, 
 // i=8 muons, i=9 neutrinos
}

void  RemSimAnalysisManager::SecondaryProtonReachingThePhantom(G4double energy)
{
 histo_proton_reaching -> fill(energy);
 //G4cout<< "proton energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryNeutronReachingThePhantom(G4double energy)
{
  histo_neutron_reaching -> fill(energy);
  //G4cout<< "neutron energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryPionReachingThePhantom(G4double energy)
{
  histo_pion_reaching -> fill(energy);
  //G4cout<< "pion energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryAlphaReachingThePhantom(G4double energy)
{
  histo_alpha_reaching -> fill(energy);
  //G4cout<< "alpha energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryPositronReachingThePhantom(G4double energy)
{
  histo_positron_reaching -> fill(energy);
  //G4cout<< "positron energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryElectronReachingThePhantom(G4double energy)
{
  histo_electron_reaching -> fill(energy);
  //G4cout<< "electron energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryGammaReachingThePhantom(G4double energy)
{ 
  histo_gamma_reaching -> fill(energy);
  //G4cout<< "gamma energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryMuonReachingThePhantom(G4double energy)
{ 
  histo_muon_reaching -> fill(energy);
  //G4cout<< "muon energy reaching the phantom: "<< energy << G4endl;
}

void  RemSimAnalysisManager::SecondaryOtherReachingThePhantom(G4double energy)
{ 
  histo_other_reaching -> fill(energy);
 //G4cout<< "other energy reaching the phantom: "<< energy << G4endl;
}

void RemSimAnalysisManager::SecondaryInVehicle(G4int i)
{
  histo_vehicle -> fill(i);
}

void RemSimAnalysisManager::finish() 
{  
  theTree -> commit();
  theTree -> close();
}
#endif











