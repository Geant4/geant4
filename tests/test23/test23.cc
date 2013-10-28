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

#include "G4Timer.hh"

// random engine/seed settings
#include "Randomize.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4GeometryManager.hh"
#include "G4StateManager.hh"

// #include "G4HadronCrossSections.hh"

#include "TstPhysListReader.hh"
#include "TstTarget.hh"

#include "TstPrimaryGeneratorAction.hh"
#include "Tst23SteppingAction.hh"

#include "Tst23Histo.hh"

// #include "Tst23ActionInit.hh"

#include "FTFP_BERT.hh"
#include "QGSP_FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "FTF_BIC.hh"

// experimental list !!!
#include "NuBeam.hh"

#include "G4ParticleTable.hh"

using namespace std;

int main(int argc, char** argv) 
{

   TstPhysListReader* theConfigReader = new TstPhysListReader();

   // Choose the Random engine
   //
   CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
   CLHEP::HepRandom::setTheSeed( 1234 ); // just something... to be on the safe side
   
   // Control on input
   if (argc < 2) 
   {
      G4cout << "Input file is not specified! Exit" << G4endl;
      exit(1);
   }
  
   // std::string input = argv[1];
   // G4STring input = argv[1];
   theConfigReader->OpenAppConfig( argv[1] );
   theConfigReader->Help();

   // Note-1: Run manager once and for all;
   //       I can NOT create/delete run manager because
   //       once it's deleated, it'll also cause deletion
   //       of several other singleton - some in the way that
   //       will prevent their re-creation in the same job;
   //       however, run manager can be somewhat re-configured
   //   
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(2);
#else
   G4RunManager* runManager = new G4RunManager();
#endif

   // Note-2: This is related to an attempt to run a loop
   //       over physics lists, changing them as we go.
   //       It turned out to be "possible with restrictions",
   //       but is very much a headache.
   //       The reason I keep this piece of code, commented,
   //       is just a written record for myself.
   //       The reason why I can NOT create/delete a list
   //       is that when a list is deleted, it EMPTIES
   //       the particle table - and the table can NOT
   //       be re-populated because instances of particles
   //       are already there in the memory, so the system
   //       will check on that and wont move any further
   //       ... thus we'll have particles everywhere in 
   //       the memory but not in the table... oops !
   //
   // std::vector<G4VModularPhysicsList*> listCollection;

   // Note-3: All other elements can be re-created in the loop,
   //         although a run with just one phys.list is best,
   //         thus there's no real need to recreate anything...
   //
   TstTarget* geom = 0;      
   TstPrimaryGeneratorAction* beam = 0;
   TstHisto* histo = 0;
   Tst23SteppingAction* stepping = 0;
   // G4VCrossSectionDataSet* cs = new G4HadronInelasticDataSet(); 
   //
   // from 4.10.b01 onwards as it require specific (new) infrastructure
   // -->Tst23ActionInit* action = 0;
        
   do 
   {

      theConfigReader->ProcessConfig();
      if ( theConfigReader->IsDone() ) break; // double protection because 
                                              // in the previous cycle it stopped at #run

      if(!G4StateManager::GetStateManager()->SetNewState(G4State_PreInit))
         G4cout << "G4StateManager PROBLEM! " << G4endl;
      
      CLHEP::HepRandom::setTheSeed( theConfigReader->GetRndmSeed() ); // reset the seed

      if ( geom ) delete geom;
      geom = new TstTarget();    
      G4ThreeVector targetSize = theConfigReader->GetTargetSize();
      geom->SetDimentions( targetSize.x(), targetSize.y(), targetSize.z() );
      geom->SetShape( theConfigReader->GetTargetShape() );      
      geom->ResetMaterial( theConfigReader->GetTargetMaterial() );
      
      runManager->SetUserInitialization( geom ); 
      runManager->GeometryHasBeenModified();
                     
// Again, this commented potion below is related to an attempt 
// to run a loop over physics lists.
// See a note earlier in the code.
/*
      if ( theConfigReader->GetPhysics() == "ftfp_bert" )
      {
	 listCollection.push_back( new FTFP_BERT() );
      } 
      else if ( theConfigReader->GetPhysics() == "qgsp_ftfp_bert" )
      {
	 listCollection.push_back( new QGSP_FTFP_BERT() );
      }
      else if ( theConfigReader->GetPhysics() == "qgsp_bert" )
      {
	 listCollection.push_back( new QGSP_BERT() );
      }
      else if (  theConfigReader->GetPhysics() == "ftf_bic" )
      {
         listCollection.push_back( new FTF_BIC() );
      } 
      else if (  theConfigReader->GetPhysics() ==  "NuBeam" )
      {
         listCollection.push_back( new TentativeNuMI() ) ;
      }
      
      int listCounter = listCollection.size() - 1;
      runManager->SetUserInitialization( listCollection[listCounter] ); 
      // runManager->PhysicsHasBeenModified();    
*/      
      
      G4VModularPhysicsList* plist = 0;
      if ( theConfigReader->GetPhysics() == "ftfp_bert" )
      {
	 plist = new FTFP_BERT();
      } 
      else if ( theConfigReader->GetPhysics() == "qgsp_ftfp_bert" )
      {
	 plist = new QGSP_FTFP_BERT();
      }
      else if ( theConfigReader->GetPhysics() == "qgsp_bert" )
      {
	 plist = new QGSP_BERT();
      }
      else if (  theConfigReader->GetPhysics() == "ftf_bic" )
      {
         plist = new FTF_BIC();
      } 
      else if (  theConfigReader->GetPhysics() ==  "NuBeam" )
      {
         plist = new NuBeam();
      }
      
      runManager->SetUserInitialization( plist );
      
      if ( beam ) delete beam;
      beam = new TstPrimaryGeneratorAction();
      beam->InitBeam( theConfigReader );
      // up to 4.9.6.p02
      runManager->SetUserAction( beam );
      
      if ( histo ) delete histo;
      histo = new Tst23Histo( theConfigReader );
      
      if ( stepping ) delete stepping;
      stepping = new Tst23SteppingAction( histo );
      // up to 4.9.6.p02
      runManager->SetUserAction( stepping );
      
/* from 4.10.b01 onwards
      if ( action ) delete action;
      action = new Tst23ActionInit();
      action->SetAct( beam );
      action->SetAct( stepping );
      
      runManager->SetUserInitialization( action );
*/
      if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))            
         G4cout << "G4StateManager PROBLEM! " << G4endl;

      runManager->InitializeGeometry();
      runManager->InitializePhysics();      
      runManager->Initialize();

      //G4double xsec = beam->GetXSecOnTarget() / millibarn;
      //G4cout << " xsec = " << xsec << G4endl;

      runManager->ConfirmBeamOnCondition();
      runManager->RunInitialization(); // this is part of BeamOn 
                                       // and needs be done (at least) to set GeomClosed status 
      
      stepping->SetTargetPtr( geom->GetTarget() );
      
      for (G4int iter=0; iter<theConfigReader->GetNEvents(); ++iter) 
      {
	 // G4cout << " event # " << iter << G4endl;
	 runManager->ProcessOneEvent( iter );
         // G4Event* event = runManager->GenerateEvent( iter );
      }
      
      histo->Write( theConfigReader->GetNEvents(), (beam->GetXSecOnTarget()/millibarn) );
      
      runManager->RunTermination();
      
   } while(!theConfigReader->IsDone());

   G4cout << "###### End of test #####" << G4endl;
      
// See notes earlier in the code, on phys.lists business
/*
   int NLists = listCollection.size();
   
   // remove all physics lists but the last one 
   // as it HAS to be removed by run manager dtor
   //
   for ( int i=0; i<NLists-1; i++ )
   {
      if ( listCollection[i] ) 
      {
         delete listCollection[i];
	 listCollection[i] = 0;
      }
   }
*/
   delete runManager; 
    
   delete theConfigReader; 
   
   return 0;

}
