//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: DMXAnalysisManager.cc
// GEANT4 tag $Name:
//
// Author: Alex Howard (a.s.howard@ic.ac.uk)
//
// History:
// -----------
// 16 Jan 2002 Alex Howard     Created
// 17 June 2002 Alex Howard    Successfully Modified to AIDA 2.2
//
// -------------------------------------------------------------------
#ifdef  G4ANALYSIS_USE

#include "DMXAnalysisManager.hh"

//#include "g4std/iomanip"

DMXAnalysisManager* DMXAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager::DMXAnalysisManager() :
  af(0), tree(0), hf(0), tpf(0), pf(0)
{
  // tree is created and booked inside book()
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager::~DMXAnalysisManager() 
{

  delete pf;
  pf=0;

 delete tpf;
  tpf=0;

  delete hf;
  hf=0;

  delete tree;
  tree=0;

  delete af;
  af = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXAnalysisManager* DMXAnalysisManager::getInstance()
{
  if (instance == 0) instance = new DMXAnalysisManager;
  return instance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void DMXAnalysisManager::book(G4String histogramfile)

{
  //  histoManager->selectStore("DMX.his");
  G4cout << " Histogramfile: " << histogramfile << G4endl;


  //build up  the  factories
  af = AIDA_createAnalysisFactory();


 //parameters for the TreeFactory
  G4bool fileExists = false;
  G4bool readOnly   = false;

  AIDA::ITreeFactory     * tf = af->createTreeFactory();

  tree = tf->create(histogramfile, "hbook", readOnly, fileExists);

  G4cout << "Tree store : " << tree->storeName() << G4endl;

  G4cout << " Booked Hbook File " << G4endl;

  //HistoFactory and TupleFactory depend on theTree
  hf = af->createHistogramFactory( *tree );
  tpf  = af->createTupleFactory(*tree );

 // ---- primary ntuple ------

  AIDA::ITuple* ntuple1 = tpf->create( "1", "Particle Source Energy", 
			     "double energy" );

  assert(ntuple1);

  // ---- secondary ntuple ------   

  AIDA::ITuple* ntuple2 = tpf->create( "2", "Scintillation Hits Info", 
				 "float Event,e_prim,tot_e,s_hits,xe_time,num_ph,avphtime,1stpart,1stparte,gamma,neutron,posi,elec,other,seed1,seed2" );

  assert(ntuple2);

  // ---- tertiary ntuple ------   

 AIDA::ITuple* ntuple3 = tpf->create( "3", "PMT Hits Info", 
				"float event, hits, xpos, ypos, zpos" );

 assert(ntuple3);

  // ---- extra ntuple ------   
  AIDA::ITuple* ntuple4 = tpf->create( "4", "Particles energy type", 
			     "float energy, NameIdx" );

  assert(ntuple4);


  // Creating an 1-dimensional histogram in the root directory of the tree

  AIDA::IHistogram1D* hEsourcep;
  hEsourcep    = hf->createHistogram1D("10","Source Energy /keV",  1000,0.,10000.);

  AIDA::IHistogram1D* hEdepp;
  hEdepp       = hf->createHistogram1D("20","Energy Deposit /keV", 1000,0.,1000.);
  
  AIDA::IHistogram1D* hEdepRecoil;
  hEdepRecoil  = hf->createHistogram1D("30","Nuclear Recoil Edep /keV", 100,0.,100.);
  
  AIDA::IHistogram1D* hNumPhLow;
  hNumPhLow    = hf->createHistogram1D("40","Number of Photons - LowE", 200,0.,200.);
  
  AIDA::IHistogram1D* hNumPhHigh;
  hNumPhHigh   = hf->createHistogram1D("50","Number of Photons - HighE", 100,0.,10000.);
  
  AIDA::IHistogram1D* hAvPhArrival;
  hAvPhArrival  = hf->createHistogram1D("60","Average Photon Arrival/ns", 200,0.,200.);
  AIDA::IHistogram1D* h1stPhArrival;
  h1stPhArrival = hf->createHistogram1D("61","1st event Photon Arrival", 200,0.,200.);
  
  AIDA::IHistogram2D* hPMTHits;
  hPMTHits    = hf->createHistogram2D("70","PMT Hit Pattern", 
			  300 ,-30.,30.,300,-30.,30.);
  AIDA::IHistogram2D* h1stPMTHit;
  h1stPMTHit  = hf->createHistogram2D("71","1st event PMT Hit Pattern", 
			     300 ,-30.,30.,300,-30.,30.);

  AIDA::IHistogram1D* hGammaEdep;
  hGammaEdep    = hf->createHistogram1D("91","Gamma Energy Deposit/keV", 1000,0.,1000.);
  AIDA::IHistogram1D* hNeutronEdep;
  hNeutronEdep  = hf->createHistogram1D("92","Neutron Ener Deposit/keV", 1000,0.,1000.);
  AIDA::IHistogram1D* hElectronEdep;
  hElectronEdep = hf->createHistogram1D("93","Electron Ener Deposit/keV",1000,0.,1000.);
  AIDA::IHistogram1D* hPositronEdep;
  hPositronEdep = hf->createHistogram1D("94","Positron Ener Deposit/keV",1000,0.,1000.);
  AIDA::IHistogram1D* hOtherEdep;
  hOtherEdep    = hf->createHistogram1D("95","Other Ener Deposit/keV", 1000,0.,1000.);
  
  delete tf;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::finish()
{

  // Committing the transaction with the tree
  G4std::cout << "Committing..." << G4std::endl;
  // write all histograms to file
  tree->commit();

  G4std::cout << "Closing the tree..." << G4std::endl;

  // close (will again commit)
  tree->close();

  // extra delete as objects are created in book() method rather than during
  // initialisation of class
  delete pf;
  delete tpf;
  delete hf;
  delete tree;
  delete af;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
    AIDA::IHistogram1D* h4 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("30") );
    h4->fill(totEnergy);  // fill(x,y,weight)     
  }

  AIDA::IHistogram1D* h2 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("40") );
  h2->fill(P_hits,10.);  // fill(x,weight) 

  AIDA::IHistogram1D* h3 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("50") );
  h3->fill(P_hits);  // fill(x,y,weight) 


  AIDA::IHistogram1D* h1 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("10") );
  h1->fill( energy_pri/keV );  // fill(x,weight)     

  AIDA::IHistogram1D* h5 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("20") );
  h5->fill( totEnergy/keV );
  
  AIDA::IHistogram1D* h6 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("60") );
  h6->fill(aveTimePmtHits/ns);  // fill(x,y,weight)     

  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("2") );

 // Fill the ntuple
  ntuple->fill( ntuple->findColumn( "Event"   ), (G4float) event_id          );
  ntuple->fill( ntuple->findColumn( "e_prim"  ), (G4float) energy_pri/keV    );
  ntuple->fill( ntuple->findColumn( "tot_e"   ), (G4float) totEnergy         );
  ntuple->fill( ntuple->findColumn( "s_hits"  ), (G4float) S_hits            );
  ntuple->fill( ntuple->findColumn( "xe_time" ), (G4float) firstLXeHitTime   );
  ntuple->fill( ntuple->findColumn( "num_ph"  ), (G4float) P_hits            );
  ntuple->fill( ntuple->findColumn( "avphtime"), (G4float) aveTimePmtHits    );
  ntuple->fill( ntuple->findColumn( "1stpart" ), (G4float) firstparticleIndex);
  ntuple->fill( ntuple->findColumn( "1stparte"), (G4float) firstParticleE    );
  ntuple->fill( ntuple->findColumn( "gamma"   ), (G4float) gamma_ev          );
  ntuple->fill( ntuple->findColumn( "neutron" ), (G4float) neutron_ev        );
  ntuple->fill( ntuple->findColumn( "posi"    ), (G4float) positron_ev       );
  ntuple->fill( ntuple->findColumn( "elec"    ), (G4float) electron_ev       );
  ntuple->fill( ntuple->findColumn( "other"   ), (G4float) other_ev          );
  ntuple->fill( ntuple->findColumn( "seed1"   ), (G4float) seed1             );
  ntuple->fill( ntuple->findColumn( "seed2"   ), (G4float) seed2             );

//ntuple->fill( ntuple->findColumn( "Event" ),  static_cast<float>(event_id) );

  //Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::analysePMTHits(G4int event, G4int i, G4double x, G4double y, G4double z)
{

  AIDA::IHistogram2D* h7 = dynamic_cast<AIDA::IHistogram2D *> ( tree->find("70") );
  h7->fill(x/mm, y/mm);  // fill(x,y,weight)     

  if (event == 0 ) {
    AIDA::IHistogram2D* h9 = dynamic_cast<AIDA::IHistogram2D *> ( tree->find("71") );
    h9->fill(x,y); // fill(x,y,weight)
  }

  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("3") );
  // Fill the secondaries ntuple
  ntuple->fill( ntuple->findColumn( "event" ), (G4float) event );
  ntuple->fill( ntuple->findColumn( "hits"  ), (G4float) i     );
  ntuple->fill( ntuple->findColumn( "xpos"  ), (G4float) x     );
  ntuple->fill( ntuple->findColumn( "ypos"  ), (G4float) y     );
  ntuple->fill( ntuple->findColumn( "zpos"  ), (G4float) z     );

  // NEW: Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow(); // check for returning true ...

  //  IHistogram1D* h3 = dynamic_cast<IHistogram1D *> ( theTree->find("30") );
  //h3->fill(ph_hits);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::analysePrimaryGenerator(G4double energy)
{

  //  IHistogram1D* h1 = dynamic_cast<IHistogram1D *> ( tree->find("10") );
  //  h1->fill(static_cast<double>(energy/keV));  // fill(x,weight) 

  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("1") );
  // Fill energy ntple:
  ntuple->fill( ntuple->findColumn( "energy" ), energy );

  // NEW: Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow(); // check for returning true ...

  // h1 and h2 are the order in which they are filled, "1" and "2"
  // are the labels associated with the histogram
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::analyseParticleSource(G4double energy, G4String name)
{

  if(name == "gamma") {
    AIDA::IHistogram1D* h11 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("91") );
    h11->fill(energy/keV);  // fill(x,weight)     
  }
  if(name == "neutron") {
    AIDA::IHistogram1D* h12 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("92") );
    h12->fill(energy/keV);  // fill(x,weight)     
  }    
    if(name == "electron") {
    AIDA::IHistogram1D* h13 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("93") );
    h13->fill(energy/keV);  // fill(x,weight)     
  }    
      if(name == "positron") {
    AIDA::IHistogram1D* h14 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("94") );
    h14->fill(energy/keV);  // fill(x,weight)     
  }    
	if(name == "other") {
    AIDA::IHistogram1D* h15 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("95") );
    h15->fill(energy/keV);  // fill(x,weight)     
  }    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::HistFirstTime(G4double time)
{   
  AIDA::IHistogram1D* h8 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("61") );
  h8->fill(time/ns);  // fill(x,y,weight)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXAnalysisManager::PlotHistos(G4bool interactive)
{   

  AIDA::IHistogram1D* h1p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("10") );
  AIDA::IHistogram1D& h1  = *h1p;  
  AIDA::IHistogram1D* h2p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("20") );
  AIDA::IHistogram1D& h2  = *h2p;  
  AIDA::IHistogram1D* h3p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("40") );
  AIDA::IHistogram1D& h3  = *h3p;  
  AIDA::IHistogram1D* h4p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("50") );
  AIDA::IHistogram1D& h4  = *h4p;  
  AIDA::IHistogram1D* h5p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("60") );
  AIDA::IHistogram1D& h5  = *h5p;  
  AIDA::IHistogram1D* h6p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("61") );
  AIDA::IHistogram1D& h6  = *h6p;  
  AIDA::IHistogram2D* h7p = dynamic_cast<AIDA::IHistogram2D *> ( tree->find("70") );
  AIDA::IHistogram2D& h7  = *h7p;  
  AIDA::IHistogram1D* h8p = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("91") );
  AIDA::IHistogram1D& h8  = *h8p;  

//   // Creating the plotter factory
//   pf = af->createPlotterFactory();
//   // Creating a plotter
//   AIDA::IPlotter* plotter = pf->create();
//   //  plotter = pf->create();

//   // Creating two regions
//   plotter->clearRegions();
//   //  plotter->createRegions(2, 2, 0); // set the current working region to the first one
//   plotter->show();

//   // Plotting the second histogram in the first region
//   plotter->plot( h1 );
//   plotter->refresh();
//   plotter->write("summary1.ps", "ps");

//   // Plotting the first histogram in the next available region
//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h2 );
//   plotter->refresh();
//   plotter->write("summary2.ps", "ps");

//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h3 );
//   plotter->refresh();
//   plotter->write("summary3.ps", "ps");

//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h4 );
//   plotter->refresh();
//   plotter->write("summary4.ps", "ps");

//   // Update the canvas on the screen
//   //  plotter->refresh();

//   //plotter->write("summary1.ps", "ps");


//   if (interactive) {
//     // Wait for the keyboard return to avoid destroying the plotter window too quickly.
//     G4cout << "Press <ENTER> to exit" << G4endl;
//     G4cin.get();
//   }

//   plotter = pf->create();
//   //  plotter->clearPage();
//   plotter->clearRegions();
//   //  plotter->createRegions(2, 2, 0); // set the current working region to the first one
//   plotter->show();
//   //  plotter->setCurrentRegion( 0 );

//   plotter->plot( h5 );
//   plotter->refresh();
//   plotter->write("summary5.ps", "ps");

//   // Plotting the first histogram in the next available region
//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h6 );
//   plotter->refresh();
//   plotter->write("summary6.ps", "ps");

//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h7 );
//   plotter->refresh();
//   plotter->write("summary7.ps", "ps");

//   //  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
//   plotter->plot( h8 );
//   plotter->refresh();
//   plotter->write("summary8.ps", "ps");

//   // Update the canvas on the screen
//   //  plotter->refresh();

//   //  plotter->write("summary2.ps", "ps");

//   if (interactive) {
//     // Wait for the keyboard return to avoid destroying the plotter window too quickly.
//     G4cout << "Press <ENTER> to exit" << G4endl;
//     G4cin.get();
//   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif







