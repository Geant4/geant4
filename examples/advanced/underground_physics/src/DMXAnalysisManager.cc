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

  ITreeFactory     * tf = af->createTreeFactory();

  tree = tf->create(histogramfile, readOnly, fileExists, "hbook");

  G4cout << "Tree store : " << tree->storeName() << G4endl;

  G4cout << " Booked Hbook File " << G4endl;

  //HistoFactory and TupleFactory depend on theTree
  hf = af->createHistogramFactory( *tree );
  tpf  = af->createTupleFactory(*tree );

 // ---- primary ntuple ------

  ITuple* ntuple1 = tpf->create( "1 Energy", "Particle Source Energy", 
			     "float energy" );

  assert(ntuple1);

  // ---- secondary ntuple ------   

  ITuple* ntuple2 = tpf->create( "2 Hits Info", "Scintillation Hits Info", 
				 "float Event,e_prim,tot_e,s_hits,xe_time,num_ph,avphtime,1stpart,1stparte,gamma,neutron,posi,elec,other,seed1,seed2" );

  assert(ntuple2);

  // ---- tertiary ntuple ------   

 ITuple* ntuple3 = tpf->create( "3 Pmt Info", "PMT Hits Info", 
				"float event, hits, xpos, ypos, zpos" );

 assert(ntuple3);

  // ---- extra ntuple ------   
  ITuple* ntuple4 = tpf->create( "4 Particle Info", "Particles energy type", 
			     "float energy, NameIdx" );

  assert(ntuple4);


  // Creating an 1-dimensional histogram in the root directory of the tree

  IHistogram1D* hEsourcep;
  hEsourcep    = hf->create1D("10","Source Energy /keV",  1000,0.,10000.);

  IHistogram1D* hEdepp;
  hEdepp       = hf->create1D("20","Energy Deposit /keV", 1000,0.,1000.);
  
  IHistogram1D* hEdepRecoil;
  hEdepRecoil  = hf->create1D("30","Nuclear Recoil Edep /keV", 100,0.,100.);
  
  IHistogram1D* hNumPhLow;
  hNumPhLow    = hf->create1D("40","Number of Photons - LowE", 200,0.,200.);
  
  IHistogram1D* hNumPhHigh;
  hNumPhHigh   = hf->create1D("50","Number of Photons - HighE", 100,0.,10000.);
  
  IHistogram1D* hAvPhArrival;
  hAvPhArrival  = hf->create1D("60","Average Photon Arrival/ns", 200,0.,200.);
  IHistogram1D* h1stPhArrival;
  h1stPhArrival = hf->create1D("61","1st event Photon Arrival", 200,0.,200.);
  
  IHistogram2D* hPMTHits;
  hPMTHits    = hf->create2D("70","PMT Hit Pattern", 
			  300 ,-30.,30.,300,-30.,30.);
  IHistogram2D* h1stPMTHit;
  h1stPMTHits = hf->create2D("71","1st event PMT Hit Pattern", 
			     300 ,-30.,30.,300,-30.,30.);

  IHistogram1D* hGammaEdep;
  hGammaEdep    = hf->create1D("91","Gamma Energy Deposit/keV", 1000,0.,1000.);
  IHistogram1D* hNeutronEdep;
  hNeutronEdep  = hf->create1D("92","Neutron Ener Deposit/keV", 1000,0.,1000.);
  IHistogram1D* hElectronEdep;
  hElectronEdep = hf->create1D("93","Electron Ener Deposit/keV",1000,0.,1000.);
  IHistogram1D* hPositronEdep;
  hPositronEdep = hf->create1D("94","Positron Ener Deposit/keV",1000,0.,1000.);
  IHistogram1D* hOtherEdep;
  hOtherEdep    = hf->create1D("95","Other Ener Deposit/keV", 1000,0.,1000.);
  
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
    IHistogram1D* h4 = dynamic_cast<IHistogram1D *> ( tree->find("30") );
    h4->fill(totEnergy);  // fill(x,y,weight)     
  }

  IHistogram1D* h2 = dynamic_cast<IHistogram1D *> ( tree->find("40") );
  h2->fill(P_hits,10.);  // fill(x,weight) 

  IHistogram1D* h3 = dynamic_cast<IHistogram1D *> ( tree->find("50") );
  h3->fill(P_hits);  // fill(x,y,weight) 


  IHistogram1D* h1 = dynamic_cast<IHistogram1D *> ( tree->find("10") );
  h1->fill( energy_pri/keV );  // fill(x,weight)     

  IHistogram1D* h5 = dynamic_cast<IHistogram1D *> ( tree->find("20") );
  h5->fill( totEnergy/keV );
  
  IHistogram1D* h6 = dynamic_cast<IHistogram1D *> ( tree->find("60") );
  h6->fill(aveTimePmtHits/ns);  // fill(x,y,weight)     

  ITuple * ntuple = dynamic_cast<ITuple *> ( tree->find("2 Hits Info") );

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

  IHistogram2D* h7 = dynamic_cast<IHistogram2D *> ( tree->find("70") );
  h7->fill(x/mm, y/mm);  // fill(x,y,weight)     

  if (event == 0 ) {
    IHistogram2D* h9 = dynamic_cast<IHistogram2D *> ( tree->find("71") );
    h9->fill(x,y); // fill(x,y,weight)
  }

  ITuple * ntuple = dynamic_cast<ITuple *> ( tree->find("3 Pmt Info") );
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

  ITuple * ntuple = dynamic_cast<ITuple *> ( tree->find("1 Energy") );
  // Fill energy ntple:
  ntuple->fill( ntuple->findColumn( "energy" ), (G4float) energy );

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
    IHistogram1D* h11 = dynamic_cast<IHistogram1D *> ( tree->find("91") );
    h11->fill(energy/keV);  // fill(x,weight)     
  }
  if(name == "neutron") {
    IHistogram1D* h12 = dynamic_cast<IHistogram1D *> ( tree->find("92") );
    h12->fill(energy/keV);  // fill(x,weight)     
  }    
    if(name == "electron") {
    IHistogram1D* h13 = dynamic_cast<IHistogram1D *> ( tree->find("93") );
    h13->fill(energy/keV);  // fill(x,weight)     
  }    
      if(name == "positron") {
    IHistogram1D* h14 = dynamic_cast<IHistogram1D *> ( tree->find("94") );
    h14->fill(energy/keV);  // fill(x,weight)     
  }    
	if(name == "other") {
    IHistogram1D* h15 = dynamic_cast<IHistogram1D *> ( tree->find("95") );
    h15->fill(energy/keV);  // fill(x,weight)     
  }    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXAnalysisManager::HistFirstTime(G4double time)
{   
  IHistogram1D* h8 = dynamic_cast<IHistogram1D *> ( tree->find("61") );
  h8->fill(time/ns);  // fill(x,y,weight)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXAnalysisManager::PlotHistos(G4bool interactive)
{   

  IHistogram1D* h1p = dynamic_cast<IHistogram1D *> ( tree->find("10") );
  IHistogram1D& h1  = *h1p;  
  IHistogram1D* h2p = dynamic_cast<IHistogram1D *> ( tree->find("20") );
  IHistogram1D& h2  = *h2p;  
  IHistogram1D* h3p = dynamic_cast<IHistogram1D *> ( tree->find("40") );
  IHistogram1D& h3  = *h3p;  
  IHistogram1D* h4p = dynamic_cast<IHistogram1D *> ( tree->find("50") );
  IHistogram1D& h4  = *h4p;  
  IHistogram1D* h5p = dynamic_cast<IHistogram1D *> ( tree->find("60") );
  IHistogram1D& h5  = *h5p;  
  IHistogram1D* h6p = dynamic_cast<IHistogram1D *> ( tree->find("61") );
  IHistogram1D& h6  = *h6p;  
  IHistogram2D* h7p = dynamic_cast<IHistogram2D *> ( tree->find("70") );
  IHistogram2D& h7  = *h7p;  
  IHistogram1D* h8p = dynamic_cast<IHistogram1D *> ( tree->find("91") );
  IHistogram1D& h8  = *h8p;  

  // Creating the plotter factory
  pf = af->createPlotterFactory();
  // Creating a plotter
  IPlotter* plotter = pf->create();
  //  plotter = pf->create();

  // Creating two regions
  plotter->clearPage();
  plotter->createRegions(2, 2, 0); // set the current working region to the first one
  plotter->show();

  // Plotting the second histogram in the first region
  plotter->plot( h1 );

  // Plotting the first histogram in the next available region
  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h2 );

  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h3 );

  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h4 );

  // Update the canvas on the screen
  plotter->refresh();

  plotter->write("summary1.ps", "ps");


  if (interactive) {
    // Wait for the keyboard return to avoid destroying the plotter window too quickly.
    G4cout << "Press <ENTER> to exit" << G4endl;
    G4cin.get();
  }

  plotter = pf->create();
  plotter->clearPage();
  plotter->createRegions(2, 2, 0); // set the current working region to the first one
  plotter->show();
  //  plotter->setCurrentRegion( 0 );

  plotter->plot( h5 );

  // Plotting the first histogram in the next available region
  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h6 );

  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h7 );

  plotter->next();    // or explicitly :  plotter->setCurrentRegion( 1 );
  plotter->plot( h8 );

  // Update the canvas on the screen
  plotter->refresh();

  plotter->write("summary2.ps", "ps");

  if (interactive) {
    // Wait for the keyboard return to avoid destroying the plotter window too quickly.
    G4cout << "Press <ENTER> to exit" << G4endl;
    G4cin.get();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif







