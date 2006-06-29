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
#include "Tst14ProcCallSA.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

#ifdef histo
#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>
#endif

Tst14ProcCallSA::Tst14ProcCallSA()
{
#ifdef histo
  // init hbook
  hbookManager = new HBookFile("proccall.his", 1);
  assert (hbookManager != 0);
#endif
}

Tst14ProcCallSA::~Tst14ProcCallSA(){
#ifdef histo
   // Write histogram file
  hbookManager->write();
#endif

  print();
  
}

void Tst14ProcCallSA::execute(const G4Step* aStep){

     
     G4Track* theTrack = aStep->GetTrack();
     if(theTrack->GetNextVolume()==0 ) return;  
     G4String particleType = theTrack->GetDefinition()->GetParticleName();
     
     G4StepPoint* postStepPoint = aStep->GetPostStepPoint(); 
     G4String procname = postStepPoint->GetProcessDefinedStep()->GetProcessName();

     G4Material* material = postStepPoint->GetMaterial();
     G4String matname = material->GetName();
     
     G4String index = procname
	              + G4String(" for ") + particleType 
	              + G4String(" in ") + matname;
     calls[index] ++; 

#ifdef histo
     G4String sIndex =  particleType + G4String(" in ") + matname;
		      
     if(start[sIndex] == 0 ) {
        G4String name = G4String("All ") + sIndex
		      + G4String(" : logE (MeV)");
       char* histoName = name.data() ;
       start[sIndex] = hbookManager->histogram(histoName,220,-6.,5.);	
     }
			
     if(hist[index] == 0 ) {
        G4String name = index + G4String(" : logE (MeV)");
        char* histoName = name.data() ;
        hist[index] = hbookManager->histogram(histoName,220,-6.,5.);
     }
     G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
     G4double stepL  = aStep->GetStepLength();
     
     if(energy>0. && stepL>1e-20) {
        
        G4double weight = 1./stepL;
        start[sIndex]->accumulate(std::log10(energy) , 1 ) ;
        hist[index]->accumulate(std::log10(energy) , weight ) ;   
     }	 
#endif     
}

void Tst14ProcCallSA::print(){

     for(intMapIter icall=calls.begin(); icall!=calls.end(); icall++){
         G4cout<<(*icall).first<<" : "<<(*icall).second<<" calls"<<G4endl;
     }

}

