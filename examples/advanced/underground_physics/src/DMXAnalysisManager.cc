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
// $Id: DMXAnalysisManager.cc
// GEANT4 tag $Name:
//
// Author: Alex Howard (alexander.howard@cern.ch)
//
// History:
// -----------
// 16 Jan 2002 Alex Howard     Created
// 17 June 2002 Alex Howard    Successfully Modified to AIDA 2.2
// 17 November 2002 Alex Howard Migrated to AIDA 3.0 and added fitting
// 22 October 2009 Luciano Pandola Changed variables names for ntuple2
//                                 (avoid problem with columns) and some 
//                                 clean-up (e.g. remove ntuple4)
//
// -------------------------------------------------------------------
#ifdef  G4ANALYSIS_USE

#include "DMXAnalysisManager.hh"

DMXAnalysisManager* DMXAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager::DMXAnalysisManager() :
  af(0), tf(0),tree(0), hf(0), tpf(0),
   ntuple1(0), ntuple2(0), ntuple3(0), 
  hEsourcep(0), hEdepp(0), hEdepRecoil(0), hNumPhLow(0), hNumPhHigh(0),
  hAvPhArrival(0), hPhArrival(0), hPMTHits(0), h1stPMTHit(0),hGammaEdep(0), 
  hNeutronEdep(0), hElectronEdep(0), hPositronEdep(0), hOtherEdep(0) 
  //,funFact(0),fitFact(0),exponFun(0),gaussFun(0),e_fitter(0),fitResult(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager::~DMXAnalysisManager() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager* DMXAnalysisManager::getInstance()
{
  if (instance == 0) instance = new DMXAnalysisManager();
  return instance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void DMXAnalysisManager::book(G4String histogramfile, G4bool)

{
  //  histoManager->selectStore("DMX.his");
  G4cout << " Histogramfile: " << histogramfile << G4endl;


  //build up  the  factories
  af = AIDA_createAnalysisFactory();


 //parameters for the TreeFactory
  G4bool fileExists = true;
  G4bool readOnly   = false;

  tf = af->createTreeFactory();

  tree = tf->create(histogramfile, "hbook", readOnly, fileExists);

  G4cout << "Tree store : " << tree->storeName() << G4endl;

  G4cout << " Booked Hbook File " << G4endl;

  //HistoFactory and TupleFactory depend on theTree
  hf = af->createHistogramFactory( *tree );
  tpf  = af->createTupleFactory(*tree );

 // ---- primary ntuple ------

  ntuple1 = tpf->create( "1", "Particle Source Energy", 
			     "double energy" );

  // ---- secondary ntuple ------   

  ntuple2 = tpf->create( "2", "Scintillation Hits Info", 
			 "float Event,e_prim,tot_e,s_hits,xe_time,num_ph,avphtime,firstpart,firstparte,gamma,neutron,posi,elec,other,seed1,seed2" );

  // ---- tertiary ntuple ------   

  ntuple3 = tpf->create( "3", "PMT Hits Info", 
				"float event, hits, xpos, ypos, zpos" );

  // Creating an 1-dimensional histogram in the root directory of the tree

  hEsourcep    = hf->createHistogram1D("10","Source Energy /keV",  1000,0.,10000.);

  hEdepp       = hf->createHistogram1D("20","Energy Deposit /keV", 1000,0.,1000.);
  
  hEdepRecoil  = hf->createHistogram1D("30","Nuclear Recoil Edep /keV", 100,0.,100.);
  
  hNumPhLow    = hf->createHistogram1D("40","Number of Photons - LowE", 200,0.,200.);
  
  hNumPhHigh   = hf->createHistogram1D("50","Number of Photons - HighE", 100,0.,10000.);
  
  hAvPhArrival  = hf->createHistogram1D("60","Average Photon Arrival/ns", 200,0.,200.);

  hPhArrival = hf->createHistogram1D("61","1st event Photon Arrival", 200,0.,200.);
  
  //2D histograms:
  hPMTHits    = hf->createHistogram2D("70","PMT Hit Pattern", 
			  300 ,-30.,30.,300,-30.,30.);

  h1stPMTHit  = hf->createHistogram2D("71","1st event PMT Hit Pattern", 
			     300 ,-30.,30.,300,-30.,30.);

  // extra 1D Histos:
  hGammaEdep    = hf->createHistogram1D("91","Gamma Energy Deposit/keV", 1000,0.,1000.);

  hNeutronEdep  = hf->createHistogram1D("92","Neutron Ener Deposit/keV", 1000,0.,1000.);

  hElectronEdep = hf->createHistogram1D("93","Electron Ener Deposit/keV",1000,0.,1000.);

  hPositronEdep = hf->createHistogram1D("94","Positron Ener Deposit/keV",1000,0.,1000.);

  hOtherEdep    = hf->createHistogram1D("95","Other Ener Deposit/keV", 1000,0.,1000.);

 //  // Creating the plotter factory
//   pf = af->createPlotterFactory();

//   if(pf && plotevent) {
//     plotter  = pf->create();
//     if(plotter) {
//       plotter->show();
//       plotter->setParameter("pageTitle","DMX Output Summary");
//     }
//     delete pf;
//   }

  // create fitting factories ..

 //  G4cout << " creating fit factories " << G4endl;

  //     funFact = af->createFunctionFactory(*tree);
//   fitFact = af->createFitFactory();

 //  G4cout << " creating Exponential Function " << G4endl;
//   exponFun = funFact->createFunctionByName("Exponential","E");
//   G4cout << " created " << G4endl;

//   G4cout << " creating Gaussian Function " << G4endl;
//   gaussFun = funFact->createFunctionByName("Gaussian","G");
//   G4cout << " created " << G4endl;

  //  assert(exponFun);

 //  //  e_fitter = fitFact->createFitter("Chi2","","printLevel=-1  errorUP=1.0");
//   G4cout << " creating fitter " << G4endl;
//   e_fitter = fitFact->createFitter();

//   //  assert(e_fitter); 

//   G4cout << " Created e_fitter " << G4endl;
  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::Init()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::Finish()
{

  // Committing the transaction with the tree
  G4cout << "Committing..." << G4endl;
  // write all histograms to file
  tree->commit();

  G4cout << "Closing the tree..." << G4endl;

  // close (will again commit)
  tree->close();

  
  delete tpf;
  delete hf;
  delete tree;
  delete tf;
  delete af;

  // extra delete as objects are created in book() method rather than during
  // initialisation of class
  //delete funFact;
  //delete fitFact;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXAnalysisManager::PulseTimeFit()
{

  // created fitter factory then fitter 

//   if(!plotter) plotter  = pf->create();
//   if(plotter) {
//     plotter->show();
//     plotter->setParameter("pageTitle","DMX Fitting");
//   }

//   G4cout << " Started Pulse Time Fit " << G4endl;
//   int i;
//   for (i = 0; i <  2 ; ++i) 
//     G4cout << "\n Fit Parameters: " << i << ":\t" 
// 	   << exponFun->parameterNames()[i] << G4endl;
//   // first fit to time histogram:
//   exponFun->setParameter("amp", 3.);
//   exponFun->setParameter("slope", -0.03);
//   G4cout << " parameters set " << G4endl;
  
//   fitResult = e_fitter->fit(*hPhArrival,*exponFun);
   
//   G4cout << " Fit Completed " << G4endl;
//   G4cout << " chi2/ndf is: " << fitResult->quality() << " / "
// 	 << fitResult->ndf() << G4endl;
  
//   // retrieve fitted parameters and errors and print them

//   for (i = 0; i <  2 ; ++i) 
//     G4cout << fitResult->fittedParameterNames()[i] << "   =    " 
// 	   << fitResult->fittedParameters()[i] << "   +/-  " 
// 	   << fitResult->errors()[i];

//   // plot function with its annotation and then histogram

//   plotter->currentRegion().plot(*exponFun,"[0,100] annotation");
//   plotter->currentRegion().plot(*hPhArrival, "overlay");
//   plotter->refresh();
//   plotter->writeToFile("ExponentialFit.ps","ps");


//   for (i = 0; i <  3 ; ++i) 
//     G4cout << "\n Fit Parameters: " << i << ":\t" 
// 	   << gaussFun->parameterNames()[i] << G4endl;
//   // first fit to time histogram:
//   gaussFun->setParameter("mean", 50.);
//   gaussFun->setParameter("sigma", 10.);
//   gaussFun->setParameter("amp", 2.);
//   G4cout << " parameters set " << G4endl;
  
//   fitResult = e_fitter->fit(*hAvPhArrival,*gaussFun);
   
//   G4cout << " Fit Completed " << G4endl;
//   G4cout << " chi2/ndf is: " << fitResult->quality() << " / "
// 	 << fitResult->ndf() << G4endl;
  
//   // retrieve fitted parameters and errors and print them

//   for (i = 0; i <  3 ; ++i) 
//     G4cout << fitResult->fittedParameterNames()[i] << "   =    " 
// 	   << fitResult->fittedParameters()[i] << "   +/-  " 
// 	   << fitResult->errors()[i];

//   // plot function with its annotation and then histogram

//   plotter->createRegions(1,1);
//   plotter->currentRegion().plot(*gaussFun,"[0,100] annotation");
//   plotter->currentRegion().plot(*hAvPhArrival, "overlay");
//   plotter->refresh();
//   plotter->writeToFile("TimeConstantFit.ps","ps");

//   G4cout << " Fit finished " << G4endl;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXAnalysisManager::analyseScintHits(G4int event_id, G4double energy_pri, G4double totEnergy, G4int S_hits, G4double firstLXeHitTime, G4int P_hits, G4double aveTimePmtHits, G4String firstparticleName, G4double firstParticleE, G4bool gamma_ev, G4bool neutron_ev, G4bool positron_ev,G4bool electron_ev,G4bool other_ev, long seed1, long seed2)

{

  G4int firstparticleIndex = 0;
  if(firstparticleName == "gamma") firstparticleIndex = 1;
  if(firstparticleName == "neutron") firstparticleIndex = 2;
  if(firstparticleName == "electron") firstparticleIndex = 3;
  if(firstparticleName == "positron") firstparticleIndex = 4;
  if(firstparticleName == "other") {
    firstparticleIndex = 5;
    hEdepRecoil->fill(totEnergy);
  }

  hNumPhLow->fill(P_hits,10.);  // fill(x,weight) 

  hNumPhHigh->fill(P_hits);  // fill(x,y,weight) 


  hEsourcep->fill( energy_pri/keV );

  hEdepp->fill( totEnergy/keV );
  
  hAvPhArrival->fill(aveTimePmtHits/ns);

 // Fill the ntuple
  ntuple2->fill( ntuple2->findColumn( "Event"   ), (G4float) event_id          );
  ntuple2->fill( ntuple2->findColumn( "e_prim"  ), (G4float) (energy_pri/keV)    );
  ntuple2->fill( ntuple2->findColumn( "tot_e"   ), (G4float) totEnergy         );
  ntuple2->fill( ntuple2->findColumn( "s_hits"  ), (G4float) S_hits            );
  ntuple2->fill( ntuple2->findColumn( "xe_time" ), (G4float) firstLXeHitTime   );
  ntuple2->fill( ntuple2->findColumn( "num_ph"  ), (G4float) P_hits            );
  ntuple2->fill( ntuple2->findColumn( "avphtime"), (G4float) aveTimePmtHits    );
  ntuple2->fill( ntuple2->findColumn( "firstpart" ), (G4float) firstparticleIndex);
  ntuple2->fill( ntuple2->findColumn( "firstparte"), (G4float) firstParticleE    );
  ntuple2->fill( ntuple2->findColumn( "gamma"   ), (G4float) gamma_ev          );
  ntuple2->fill( ntuple2->findColumn( "neutron" ), (G4float) neutron_ev        );
  ntuple2->fill( ntuple2->findColumn( "posi"    ), (G4float) positron_ev       );
  ntuple2->fill( ntuple2->findColumn( "elec"    ), (G4float) electron_ev       );
  ntuple2->fill( ntuple2->findColumn( "other"   ), (G4float) other_ev          );
  ntuple2->fill( ntuple2->findColumn( "seed1"   ), (G4float) seed1             );
  ntuple2->fill( ntuple2->findColumn( "seed2"   ), (G4float) seed2             );

  //Values of attributes are prepared; store them to the nTuple:
  ntuple2->addRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::analysePMTHits(G4int event, G4int i, G4double x, G4double y, G4double z)
{

  hPMTHits->fill(x/mm, y/mm);  // fill(x,y,weight)     

  if (event == 0 ) {
    h1stPMTHit->fill(x,y);
  }

  // Fill the ntuple
  ntuple3->fill( ntuple3->findColumn( "event" ), (G4float) event );
  ntuple3->fill( ntuple3->findColumn( "hits"  ), (G4float) i     );
  ntuple3->fill( ntuple3->findColumn( "xpos"  ), (G4float) x     );
  ntuple3->fill( ntuple3->findColumn( "ypos"  ), (G4float) y     );
  ntuple3->fill( ntuple3->findColumn( "zpos"  ), (G4float) z     );

  // NEW: Values of attributes are prepared; store them to the nTuple:
  ntuple3->addRow(); 

}

void DMXAnalysisManager::analysePrimaryGenerator(G4double energy)
{

  // Fill energy ntple:
  ntuple1->fill( ntuple1->findColumn( "energy" ), energy );

  // NEW: Values of attributes are prepared; store them to the nTuple:
  ntuple1->addRow(); 

}

void DMXAnalysisManager::analyseParticleSource(G4double energy, G4String name)
{

  if(name == "gamma") {
    hGammaEdep->fill(energy/keV);
  }
  if(name == "neutron") {
    hNeutronEdep->fill(energy/keV);  // fill(x,weight)     
  }    
    if(name == "electron") {
      hElectronEdep->fill(energy/keV);  // fill(x,weight)     
  }    
    if(name == "positron") {
      hPositronEdep->fill(energy/keV);  // fill(x,weight)     
  }    
    if(name == "other") {
      hOtherEdep->fill(energy/keV);  // fill(x,weight)     
  }    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::HistTime(G4double time)
{   
  hPhArrival->fill(time/ns);  // fill(x,y,weight)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::PlotHistosInter(G4int) 
{
  // // The plotter is updated only if there are some hits in the event
//   if(!flag) return;
//   //  Set the plotter ; set the number of regions and attach histograms
//   // to plot for each region.
//   //  It is done here, since then EndOfRun set regions
//   // for paper output.
//   if(plotter) {
//       plotter->currentRegion().plot(*hNumPhLow);
//       plotter->refresh();
//   }
}

void DMXAnalysisManager::PlotHistos(G4bool)
{   
  /*
  if(!plotter) plotter  = pf->create();
  plotter->show();
  plotter->setParameter("pageTitle","DMX Output Summary");

  //  G4int n = 1;
  if(plotter) {
    // We set one single region for the plotter
    // We now print the histograms, each one in a separate file
    plotter->createRegions(1,1);
    plotter->currentRegion().plot(*hNumPhHigh);
    plotter->refresh();
    plotter->writeToFile("NumberPhotons.ps","ps");
    
    plotter->createRegions(1,1);
    plotter->currentRegion().plot(*hEdepp);
    plotter->refresh();
    plotter->writeToFile("EnergyDeposit.ps","ps");

    plotter->createRegions(1,1);
    plotter->currentRegion().plot(*hEsourcep);
    plotter->refresh();
    plotter->writeToFile("SourceEnergy.ps","ps");
    
    plotter->createRegions(1,1);      
    plotter->currentRegion().plot(*hNumPhLow);
    plotter->refresh();
    plotter->writeToFile("NumberPhotonsLow.ps","ps");
  }
    //   // Creating two regions
    //   plotter->clearRegions();
    //   //  plotter->createRegions(2, 2, 0); // set the current working region to the first one
    //   plotter->show();
    
    if (interactive) {
      // Wait for the keyboard return to avoid destroying the plotter window too quickly.
      G4cout << "Press <ENTER> to exit" << G4endl;
      G4cin.get();
    }
  */
}

#endif
