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

