// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "Tst20ProcCallSA.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

Tst20ProcCallSA::Tst20ProcCallSA()
{ }

Tst20ProcCallSA::~Tst20ProcCallSA(){
  print();
  
}

void Tst20ProcCallSA::execute(const G4Step* aStep){

     
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

}

void Tst20ProcCallSA::print(){

     for(intMapIter icall=calls.begin(); icall!=calls.end(); icall++){
         G4cout<<(*icall).first<<" : "<<(*icall).second<<" calls"<<G4endl;
     }

}

