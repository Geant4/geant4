// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "Tst14ProcCallSA.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

Tst14ProcCallSA::Tst14ProcCallSA():
nparticles(0), 
nprocesses(0),
nmaterials(0)
{}

Tst14ProcCallSA::~Tst14ProcCallSA(){

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
     
     if(particles[particleType] == 0) particles[particleType] = ++nparticles;
     if(processes[procname]     == 0) processes[procname] = ++nprocesses;
     if(materials[matname]      == 0) materials[matname] = ++nmaterials;
     
     G4int index = 100*materials[matname]
                   +10*particles[particleType] 
		   +   processes[procname];
     calls[index] ++; 
}

void Tst14ProcCallSA::print(){

     for(iter ipart=particles.begin(); ipart!=particles.end(); ipart++){
         G4cout<<G4endl<<"Particle "<<(*ipart).first<<G4endl;

	 G4cout<<G4std::setw(30)<<" ";
	 for(iter imat=materials.begin(); imat!=materials.end(); imat++){
	     	 G4cout<<G4std::setw(10)<<(*imat).first;
	 }
	 G4cout<<G4endl;
	 
         for(iter iproc=processes.begin(); iproc!=processes.end(); iproc++){
	     G4cout<<G4std::setw(20) <<(*iproc).first
	           <<" called ";
	           for(iter imat=materials.begin(); imat!=materials.end(); imat++){
	               G4int index= 100*(*imat).second
		                    +10*(*ipart).second 
				    +   (*iproc).second ;
	     	       G4cout<<G4std::setw(10)<<calls[index];
	           }
	     G4cout<<" times"<<G4endl;
	 }
     }
}
