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
// $Id: G4LEInelasticTest.cc,v 1.5 2002-07-10 07:53:19 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"
 
#include "G4ProtonInelasticProcess.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4VLongitudinalStringDecay.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSMFragmentation.hh"

#include "G4DynamicParticle.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"

#include "G4ExcitationHandler.hh"
 
 // forward declarations
 
 G4int sortEnergies( const double Px, const double Py, const double Pz,
                     const double Ekin, double* sortedPx, double* sortedPy,
                     double* sortedPz, double* sortedE );
 
 // here comes the code
 
 G4int sortEnergies( const double Px, const double Py, const double Pz,
                     const double Ekin, double* sortedPx, double* sortedPy,
                     double* sortedPz, double* sortedE)
  {
    for( int i=0; i<10; ++i ) {
      if( abs(Ekin) > sortedE[i] ) {
        sortedE[i]  = Ekin;
        sortedPx[i] = Px;
        sortedPy[i] = Py;
        sortedPz[i] = Pz;
        return 1;
      }
    }
    return 0;
  }
 
 int main()
  {
    G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4std::ofstream outFile( "InInelasticAlpha.listing.GetMeanFreePath", G4std::ios::out);
    outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4std::ofstream outFile1( "InInelasticAlpha.listing.DoIt", G4std::ios::out);
    outFile1.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl, nIso;
    
 // constructing the particles
 
 G4LeptonConstructor aC1;
 G4BaryonConstructor aC2;
 G4MesonConstructor aC3;
 G4IonConstructor aC4;
 
 aC1.ConstructParticle();
 aC2.ConstructParticle();
 aC3.ConstructParticle();
 aC4.ConstructParticle();

    G4int numberOfMaterials=1;
    G4Material* theMaterials[2000];    
    
     G4cout << "Please specify A, Z"<<G4endl;
     G4cin >> a >> iz;
     G4Material *theMaterial = new G4Material(name="stuff", density=18.95*g/cm3, nEl=1);
     G4Element  *theElement  = new G4Element (name="stuff", symbol="Z", iz, a*g/mole);
     theMaterial->AddElement( theElement, 1 );
     theMaterials[25] = theMaterial;

    G4cout << "Please enter material number"<<G4endl;
    G4int inputNumber = 25;
    theMaterials[0]=theMaterials[inputNumber];
    G4cout << "Active Material = " << theMaterials[0]->GetName()<<G4endl;
    G4cin >> inputNumber;
    // ----------- here all material have been defined -----------
    
//    G4Element::DumpInfo(); 
//    G4Material::DumpInfo();
    
    // ----------- the following is needed for building a track...... ------------
    
    static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    G4int imat = 0;   
    G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);
    
    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                        (*theMaterialTable)[imat],
                                                        "LFrame", 0, 0, 0);
    
    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                                     "PFrame",LogicalFrame,0,false,0);
    G4RotationMatrix theNull;
    G4ThreeVector theCenter(0,0,0);
    G4GRSVolume * theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
    G4TouchableHandle theTouchableHandle(theTouchable);
    // ----------- now get all particles of interest ---------
   G4int numberOfParticles = 1;
   G4ParticleDefinition* theParticles[1];
   G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
   theParticles[0]=theProton;
   
   //------ here all the particles are Done ----------
   G4cout << "Done with all the particles" << G4endl;
   G4cout << "Starting process definitions" << G4endl;
   //--------- Processes definitions ---------
   G4ProtonInelasticProcess* theProcesses[1];
      
   G4ProcessManager* theProtonProcessManager = new G4ProcessManager(theProton);
   theProton->SetProcessManager(theProtonProcessManager);
   G4ProtonInelasticProcess theInelasticProcess; 
//    G4ExcitationHandler theHandler;
//   G4PreCompoundModel theProtonModel(&theHandler);
    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
    G4QGSModel<G4QGSParticipants> * theStringModel = new G4QGSModel<G4QGSParticipants>;
    G4StringChipsParticleLevelInterface * theCascade = new G4StringChipsParticleLevelInterface;
    theTheoModel->SetTransport(theCascade);
    theTheoModel->SetHighEnergyGenerator(theStringModel);
    theTheoModel->SetMinEnergy(20*GeV);
    theTheoModel->SetMaxEnergy(100*TeV);
    G4VLongitudinalStringDecay * theFragmentation;
    theFragmentation = new G4QGSMFragmentation;
    G4ExcitedStringDecay theStringDecay(theFragmentation);
    theStringModel->SetFragmentationModel(&theStringDecay);


   theInelasticProcess.RegisterMe(theTheoModel);
//   theInelasticProcess.RegisterMe(&theProtonModel);
   theProtonProcessManager->AddDiscreteProcess(&theInelasticProcess);
   theProcesses[0] = &theInelasticProcess;
   
   G4ForceCondition* condition = new G4ForceCondition;
   *condition = NotForced;

   G4cout << "Done with all the process definitions"<<G4endl;
   //   G4cout << "Building the CrossSection Tables. This will take a while"<<G4endl;
   
   // ----------- define energies of interest ----------------
      
   // ----------- needed to Build a Step, and a track -------------------
   
   G4ParticleMomentum theDirection(1.,0.,0.);
   G4ThreeVector aPosition(0.,0.,0.);
   G4double aTime = 0.0;
   
   // --------- Test the GetMeanFreePath
   
   G4StepPoint aStepPoint;
   G4Step aStep;
   aStep.SetPreStepPoint(&aStepPoint);
   G4double meanFreePath;
   G4double incomingEnergy;
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
   G4cout << "Please enter the Proton energy"<<G4endl;
   G4cin >> incomingEnergy;
   G4int errorOne;
   G4cout << "Please enter the problematic event number"<<G4endl;
   G4cin >> errorOne;
   for (i=0; i<numberOfParticles; i++)
   {
     outFile << G4endl
             << "New particle type: " << theParticles[i]->GetParticleName()
             << " " << i << G4endl;
     for ( G4int k=0; k<numberOfMaterials; k++)
     {
       outFile << G4endl << "Entering Material " << theMaterials[0]->GetName()
               << " for particle " << theParticles[i]->GetParticleName() << G4endl;
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
           aTrack->SetTouchableHandle(theTouchableHandle);
           aStep.SetTrack( aTrack );
           aStepPoint.SetTouchableHandle(theTouchableHandle);
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
	   if(hpw == errorOne)
	   {
	      hpw --;
	      hpw ++;
	   }
           G4cout << "Last chance before DoIt call: "
//                << theProtonModel.GetNiso()
                <<G4endl;
           aFinalState = (G4ParticleChange*)  (theProcesses[i]->PostStepDoIt( *aTrack, aStep ));
           G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries();
           G4double theFSEnergy = aFinalState->GetEnergyChange();
           G4ThreeVector * theFSMomentum= const_cast<G4ThreeVector *>
	                                  (aFinalState->GetMomentumChange());
           G4cout << "FINAL STATE = "<<theFSEnergy<<" ";
           G4cout <<*theFSMomentum<<G4endl;
           G4Track * second;
           G4DynamicParticle * aSec;
           G4int isec;
	   G4double QValue = 0;
	   G4double QValueM1 = 0;
	   G4double QValueM2 = 0;
	   if(aFinalState->GetStatusChange() == fAlive)
	   {
	     QValue += aFinalState->GetEnergyChange();
	     if(theParticles[i]->GetBaryonNumber()<0) QValue+= 2.*theParticles[i]->GetPDGMass();
	   }
           for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
           {
             second = aFinalState->GetSecondary(isec);
             aSec = const_cast<G4DynamicParticle *> (second->GetDynamicParticle());
             G4cout << "SECONDARIES info";
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
             G4cout << aSec->GetDefinition()->GetPDGEncoding()<<" ";
	     G4cout << (1-isec)*aFinalState->GetNumberOfSecondaries();
	     G4cout << G4endl;
	     if(aSec->GetDefinition()->GetBaryonNumber()>0)
	     {
               QValue += aSec->GetKineticEnergy();
	       if(isec!=0) QValueM1 += aSec->GetKineticEnergy();
	       if(isec>1) QValueM2 += aSec->GetKineticEnergy();
	     }
	     else if(aSec->GetDefinition()->GetBaryonNumber() == 0)
	     {
               QValue += aSec->GetTotalEnergy();
	       if(isec!=0) QValueM1 += aSec->GetTotalEnergy();
	       if(isec>1) QValueM2 += aSec->GetKineticEnergy();
	     }
	     else
	     {
               QValue += aSec->GetTotalEnergy();
               QValue += aSec->GetDefinition()->GetPDGMass();
	       if(isec!=0) QValueM1 += aSec->GetTotalEnergy();
	       if(isec!=0) QValueM1 += aSec->GetDefinition()->GetPDGMass();
	       if(isec>1) QValueM2 += aSec->GetTotalEnergy();
	       if(isec>1) QValueM2 += aSec->GetDefinition()->GetPDGMass();
	       G4cout << "found an anti !!!" <<G4endl;
	     }
	     delete second;
           }
	   G4cout << "QVALUE = "<<QValue<<" "<<QValueM1<<" "<<QValueM2<<G4endl;
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   return EXIT_SUCCESS;
}


