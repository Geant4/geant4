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

#include "TstDiscreteProcessReader.hh"

#include "Tst75ExecProcessLevel.hh"

#include "G4VParticleChange.hh"
#include "G4ParticleTable.hh"
#include "G4HadronCrossSections.hh"

#include "Tst75Histo.hh"

#include "TStopwatch.h"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

int main(int argc, char** argv) 
{

   TstDiscreteProcessReader* theConfigReader = new TstDiscreteProcessReader();
   
   // Input control
   //
   if (argc < 2) 
   {
      G4cout << "Input file is not specified! Exit" << G4endl;
      exit(1);
   }
  
   theConfigReader->OpenAppConfig( argv[1] );
   theConfigReader->Help();
  
   do 
   {

      theConfigReader->ProcessConfig();
      if ( theConfigReader->IsDone() ) break; // double protection because 
                                              // in the previous cycle it stopped at #run
           
      ExecProcessLevel* exec = new Tst75ExecProcessLevel( theConfigReader );
      
      TstHisto* histo = 0;
      histo = new Tst75Histo( theConfigReader );
      if ( !histo )
      {
         G4cout << " Histograms NOT booked" << G4endl;
         exit(1);
      }
      
      G4VParticleChange* aChange = 0;
      
      G4int NEvts = theConfigReader->GetNEvents();
      
      TStopwatch timer;
      timer.Start();
      
      for (G4int iter=0; iter<NEvts; ++iter) 
      {

         aChange = exec->DoEvent();
	 	 
	 histo->FillEvt( aChange, exec->GetBeam()->GetLabV(), exec->GetBeam()->GetLabP() );
	 
         G4int nsec = aChange->GetNumberOfSecondaries();
         for (G4int i=0; i<nsec; ++i) 
         {   
            delete aChange->GetSecondary(i);
         } 
         aChange->Clear();

      } // end loop over events
      
      timer.Stop();
      G4cout << " CPU = " << timer.CpuTime() << G4endl;
      G4cout << " Real Time = " << timer.RealTime() << G4endl;
      
      histo->Write( theConfigReader->GetNEvents(), (exec->GetXSecOnTarget()/microbarn) );
      
      if ( histo ) delete histo; // delete histograms after writing them out
      
      delete exec;

   } while(!theConfigReader->IsDone()); // end loop over "runs"

   G4cout << "###### End of test #####" << G4endl;
  
   delete theConfigReader; 
   
   return 0;

}
