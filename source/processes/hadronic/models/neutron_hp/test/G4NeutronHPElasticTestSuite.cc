// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPElasticTestSuite.cc,v 1.2 1999-12-15 14:53:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4HadronElasticProcess.hh"
 
#include "G4LElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPIsoData.hh"

#include "G4DynamicParticle.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

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
    G4std::ofstream outFile( "InelasticAlpha.listing.GetMeanFreePath", G4std::ios::out);
    outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4std::ofstream outFile1( "InelasticAlpha.listing.DoIt", G4std::ios::out);
    outFile1.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
    G4int numberOfMaterials=1;
    G4Material* theMaterials[23];
    
     G4Material *theGe = new G4Material(name="Ge72", density=18.95*g/cm3, nEl=1);
      G4Element *elGe  = new G4Element(name="Ge72", symbol="Ge72", iz=32, a=72.*g/mole);
      theGe->AddElement( elGe, 1 );
     theMaterials[1] = theGe;
    
     G4Material *theTin = new G4Material(name="Tin126", density=18.95*g/cm3, nEl=1);
      G4Element *elTin  = new G4Element(name="Tin126", symbol="Tin126", iz=50, a=126.*g/mole);
      theTin->AddElement( elTin, 1 );
     theMaterials[2] = theTin;
    
     G4Material *theHe = new G4Material(name="He4", density=18.95*g/cm3, nEl=1);
      G4Element *elHe  = new G4Element(name="He4", symbol="He4", iz=2, a=4.*g/mole);
      theHe->AddElement( elHe, 1 );
     theMaterials[3] = theHe;
    
     G4Material *theU = new G4Material(name="U238", density=18.95*g/cm3, nEl=1);
      G4Element *elU  = new G4Element(name="U238", symbol="U238", iz=92, a=238.*g/mole);
      theU->AddElement( elU, 1 );
     theMaterials[4] = theU;
    
     G4Material *theAl = new G4Material(name="Al27", density=18.95*g/cm3, nEl=1);
      G4Element *elAl  = new G4Element(name="Al27", symbol="Al27", iz=13, a=27.*g/mole);
      theAl->AddElement( elAl, 1 );
     theMaterials[5] = theAl;
    
    G4cout << "Please enter material number"<<G4endl;
    G4int inputNumber;
    G4cin >> inputNumber;
    theMaterials[0]=theMaterials[inputNumber];
    G4cout << "Active Material = " << theMaterials[0]->GetName()<<G4endl;
    G4cin >> inputNumber;
    // ----------- here all material have been defined -----------
    
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
   G4ParticleDefinition* theNeutron = G4Neutron::NeutronDefinition();
   theParticles[0]=theNeutron;
   
   //------ here all the particles are Done ----------
   G4cout << "Done with all the particles" << G4endl;
   G4cout << "Starting process definitions" << G4endl;
   //--------- Processes definitions ---------
   G4HadronElasticProcess* theProcesses[1];
      
   
   G4ProcessManager* theNeutronProcessManager = new G4ProcessManager(theNeutron);
   theNeutron->SetProcessManager(theNeutronProcessManager);
   G4HadronElasticProcess theElasticProcess; 
   G4LElastic theElastic;
   G4NeutronHPElastic theNeutronHPElastic;
   G4cout << "Elastic instanciated!!!"<<G4endl;
   theElasticProcess.RegisterMe(&theElastic);
   theElastic.SetMinEnergy(21*MeV);
   theElasticProcess.RegisterMe(&theNeutronHPElastic);
   theNeutronProcessManager->AddDiscreteProcess(&theElasticProcess);
   theProcesses[0] = &theElasticProcess;
   
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
   
   G4Step aStep;
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
   G4cout << "Please enter the neutron energy"<<G4endl;
   G4cin >> incomingEnergy;
   for (i=0; i<numberOfParticles; i++)
   {
     outFile << G4endl
             << "New particle type: " << theParticles[i]->GetParticleName()
             << " " << i << G4endl;
     for ( G4int k=0; k<numberOfMaterials; k++)
     {
       outFile << G4endl << "Entering Material " << theMaterials[k]->GetName()
               << " for particle " << theParticles[i]->GetParticleName() << G4endl;
       LogicalFrame->SetMaterial(theMaterials[k]); 
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
//                << theNeutronHPElastic.GetNiso()
                <<G4endl;
           aFinalState = (G4ParticleChange*)  (theProcesses[i]->PostStepDoIt( *aTrack, aStep ));
           G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries();
           G4double theFSEnergy = aFinalState->GetEnergyChange();
           G4ThreeVector * theFSMomentum= aFinalState->GetMomentumChange();
           G4cout << "FINAL STATE = "<<theFSEnergy<<" ";
           G4cout <<*theFSMomentum<<G4endl;
           G4Track * second;
           G4DynamicParticle * aSec;
           G4cout << "SECONDARIES info";
           G4int isec;
           for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
           {
             second = aFinalState->GetSecondary(isec);
             aSec = second->GetDynamicParticle();
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
             delete second;
           }
           G4cout << G4endl;
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   return EXIT_SUCCESS;
}


