// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: CHIPS_p_pBar_AtRest_Integration.cc,v 1.1 2000-08-19 12:14:45 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4ProtonAntiProtonAtRestChips.hh"
 
#include "G4DynamicParticle.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
 
 int main()
  {
    G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
 // constructing the particles
 
 G4LeptonConstructor aC1;
 G4BaryonConstructor aC2;
 G4MesonConstructor aC3;
 G4IonConstructor aC4;
 G4BosonConstructor aC5;
 
 aC1.ConstructParticle();
 aC2.ConstructParticle();
 aC3.ConstructParticle();
 aC4.ConstructParticle();
 aC5.ConstructParticle();

    G4int numberOfMaterials=1;
    G4Material* theMaterials[2000];
    
        G4Material *theH = new G4Material(name="Hydrogen", density=1.032*g/cm3, nEl=1);
        G4Element *elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
        theH->AddElement( elH, 1 );
        theMaterials[1] = theH;
    
    G4cout << "Please enter material number"<<G4endl;
    G4int inputNumber;
    G4cin >> inputNumber;
    theMaterials[0]=theMaterials[inputNumber];
    G4cout << "Active Material = " << theMaterials[0]->GetName()<<G4endl;
    G4cin >> inputNumber;
    
    // ----------- the following is needed for building a track...... ------------
    
    static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    G4int imat = 0;   
    G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);
    
    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                        (*theMaterialTable)(imat),
                                                        "LFrame", 0, 0, 0);
    
    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                                     "PFrame",LogicalFrame,0,false,0);
    G4RotationMatrix theNull;
    G4ThreeVector theCenter(0,0,0);
    G4GRSVolume * theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
    // ----------- now get all particles of interest ---------
   G4int numberOfParticles = 1;
   G4ParticleDefinition* theParticles[1];
   G4ParticleDefinition* theAntiProton = G4AntiProton::AntiProtonDefinition();
   theParticles[0]=theAntiProton;
   
   //------ here all the particles are Done ----------
   G4cout << "Done with all the particles" << G4endl;
   G4cout << "Starting process definitions" << G4endl;
   //--------- Processes definitions ---------
   G4VRestProcess* theProcesses[1];
      
   G4ProcessManager* theAntiProtonProcessManager = new G4ProcessManager(theAntiProton);
   theAntiProton->SetProcessManager(theAntiProtonProcessManager);
   G4ProtonAntiProtonAtRestChips theProcess; 
   G4cout << "Inelastic instanciated!!!"<<G4endl;
   theAntiProtonProcessManager->AddDiscreteProcess(&theProcess);
   theProcesses[0] = &theProcess;
   
   G4ForceCondition* condition = new G4ForceCondition;
   *condition = NotForced;

   G4cout << "Done with all the process definitions"<<G4endl;
   //   G4cout << "Building the CrossSection Tables. This will take a while"<<G4endl;
   
   // ----------- define energies of interest ----------------
   
   int numberOfEnergies = 5;
   G4double theEnergies[5] = { 0.5*MeV, 1.0*MeV, 3.0*MeV, 12.0*MeV, 20.0*MeV };
   
   // ----------- needed to Build a Step, and a track -------------------
   
   G4ParticleMomentum theDirection(0.,0.,1.);
   G4ThreeVector aPosition(0.,0.,0.);
   G4double aTime = 0.0;
   
   // --------- Test the GetMeanFreePath
   
   G4StepPoint aStepPoint;
   G4Step aStep;
   aStep.SetPreStepPoint(&aStepPoint);
   G4double meanFreePath;
   G4double incomingEnergy = 0;
   G4int k, i, l, hpw = 0;   
   // --------- Test the PostStepDoIt now, 10 events each --------------
   G4cout << "Entering the DoIt loops!!!!!"<< G4endl;
   G4ParticleChange* aFinalState;
   G4cout << "Test the DoIt: please enter the number of events"<<G4endl;
   G4int ll0;
   G4cin >> ll0;
//   G4cout <<"Now debug the DoIt: enter the problem event number"<< G4endl;
   G4int debugThisOne=1;
//   G4cin >> debugThisOne;
   for (i=0; i<numberOfParticles; i++)
   {
     for ( G4int k=0; k<numberOfMaterials; k++)
     {
       LogicalFrame->SetMaterial(theMaterials[0]); 
//       for( G4int j=0; j<numberOfEnergies; ++j )
int j = 0;
       {
         for( G4int l=0; l<ll0; ++l )
         {
//           incomingEnergy=theEnergies[j];
           G4DynamicParticle* aParticle =
             new G4DynamicParticle( theParticles[i], theDirection, incomingEnergy );
           G4Track* aTrack = new G4Track( aParticle, aTime, aPosition );
           aTrack->SetTouchable(theTouchable);
	   aStep.SetTrack( aTrack );
           aStepPoint.SetTouchable(theTouchable);
	   aStepPoint.SetMaterial(theMaterials[0]);
           aStep.SetPreStepPoint(&aStepPoint);
	   aStep.SetPostStepPoint(&aStepPoint);
	   aTrack->SetStep(&aStep);
           ++hpw;
           if(hpw == 1000*(hpw/1000))
           G4cerr << "FINAL EVENTCOUNTER=" << hpw
                << " current energy: " << incomingEnergy
                << " of particle " << aParticle->GetDefinition()->GetParticleName() 
                << " in material " << theMaterials[k]->GetName() << G4endl;
           if (hpw==debugThisOne)
           {
            debugThisOne+=0;
           }
           G4cout << "Last chance before DoIt call: "
//                << thePBarInelastic.GetNiso()
                <<G4endl;
           aFinalState = (G4ParticleChange*)  (theProcesses[i]->AtRestDoIt( *aTrack, aStep ));
           G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries();
           G4double theFSEnergy = aFinalState->GetEnergyChange();
           const G4ThreeVector * theFSMomentum= aFinalState->GetMomentumChange();
           G4cout << "FINAL STATE = "<<theFSEnergy<<" ";
           G4cout <<*theFSMomentum<<G4endl;
           G4Track * second;
           G4DynamicParticle * aSec;
           G4int isec;
           for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
           {
             second = aFinalState->GetSecondary(isec);
             aSec = second->GetDynamicParticle();
             G4cout << "SECONDARIES info";
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
	     G4cout << (1-isec)*aFinalState->GetNumberOfSecondaries();
	     G4cout << G4endl;
             delete second;
           }
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   G4cout << "CHIPS_p_pBar_Integration terminated successfully"<<endl;
   return EXIT_SUCCESS;
}


