// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPInelasticTestSuite.cc,v 1.1 1999-01-08 16:33:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#define G4DEBUG
#include "../src/G4NeutronHPEnAngCorrelation.cc"
#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>
 
#include "../src/G4NeutronHPIsoData.cc"
#include "G4Material.hh"
#include "G4Isotope.hh"

#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4NeutronInelasticProcess.hh"
 
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPIsoData.hh"

#include "G4DynamicParticle.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

#include "g4templates.hh"
#include "G4NeutronHPChannel.hh"
#include "G4Neutron.hh"

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
    G4cout.setf( ios::scientific, ios::floatfield );
    ofstream outFile( "InelasticAlpha.listing.GetMeanFreePath", ios::out);
    outFile.setf( ios::scientific, ios::floatfield );
    ofstream outFile1( "InelasticAlpha.listing.DoIt", ios::out);
    outFile1.setf( ios::scientific, ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;

    G4int numberOfMaterials=1;
    G4Material* theMaterials[23];

    cout << "Please enter Isotope A and Z"<<endl;
    cout << "A:"<<endl;
    G4int theA;
    cin >> theA;
    cout << "Z:"<<endl;
    G4int theZ;
    cin >> theZ;
    G4Material *theMat = new G4Material(name="Mat", density=18.95*g/cm3, nEl=1);
    G4Element *elMat  = new G4Element(name="Mat", symbol="M", 1);
    G4Isotope* isoMat = new G4Isotope(name="U238", theZ, theA, a=238.03*g/mole);
    elMat->AddIsotope(isoMat, 100.*perCent);
    theMat->AddElement( elMat, 1 );
    theMaterials[0]=theMat;

    G4int inputNumber;
    G4cout << "Active Material = " << theMaterials[0]->GetName()<<endl;
    cin >> inputNumber;
    
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
   G4cout << "Done with all the particles" << endl;
   G4cout << "Starting process definitions" << endl;
   //--------- Processes definitions ---------
   G4NeutronInelasticProcess* theProcesses[1];
      
   G4ProcessManager* theNeutronProcessManager = new G4ProcessManager(theNeutron);
   theNeutron->SetProcessManager(theNeutronProcessManager);
   G4NeutronInelasticProcess theInelasticProcess; 
   G4NeutronHPInelastic theNeutronHPInelastic;
   theNeutronHPInelastic.SetMaxEnergy(40*MeV);
   G4cout << "Inelastic instanciated!!!"<<endl;
   theInelasticProcess.RegisterMe(&theNeutronHPInelastic);
   theNeutronProcessManager->AddDiscreteProcess(&theInelasticProcess);
   theProcesses[0] = &theInelasticProcess;
   
   G4ForceCondition* condition = new G4ForceCondition;
   *condition = NotForced;

   G4cout << "Done with all the process definitions"<<endl;
   //   G4cout << "Building the CrossSection Tables. This will take a while"<<endl;
   
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
   G4cout << "Entering the DoIt loops!!!!!"<< endl;
   G4ParticleChange* aFinalState;
   G4cout << "Test the DoIt: please enter the number of events"<<endl;
   G4int ll0;
   cin >> ll0;
   G4cout <<"Now debug the DoIt: enter the problem event number"<< endl;
   G4int debugThisOne=1;
   cin >> debugThisOne;
   G4cout << "Please enter the neutron energy"<<endl;
   cin >> incomingEnergy;
   for (i=0; i<numberOfParticles; i++)
   {
     outFile << endl
             << "New particle type: " << theParticles[i]->GetParticleName()
             << " " << i << endl;
     for ( G4int k=0; k<numberOfMaterials; k++)
     {
       outFile << endl << "Entering Material " << theMaterials[k]->GetName()
               << " for particle " << theParticles[i]->GetParticleName() << endl;
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
           if(hpw==1000*(hpw/1000))
           G4cerr << "FINAL EVENTCOUNTER=" << hpw
                << " current energy: " << incomingEnergy
                << " of particle " << aParticle->GetDefinition()->GetParticleName() 
                << " in material " << theMaterials[k]->GetName() << endl;
           if (hpw==debugThisOne)
           {
             debugThisOne+=0;
           }
           G4cout << "Last chance before DoIt call: "
//                << theNeutronHPInelastic.GetNiso()
                <<endl;
           aFinalState =
	   (G4ParticleChange*) (theProcesses[i]->PostStepDoIt( *aTrack, aStep ));
           G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries();
           G4double theFSEnergy = aFinalState->GetEnergyChange();
           G4ThreeVector * theFSMomentum= aFinalState->GetMomentumChange();
           G4cout << "FINAL STATE = "<<theFSEnergy<<" ";
           G4cout <<*theFSMomentum<<endl;
           G4Track * second;
           G4DynamicParticle * aSec;
           G4int isec, nNeutrons;
           nNeutrons = 0;
           for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
           {
             G4cout << "SECONDARIES info";
             second = aFinalState->GetSecondary(isec);
             aSec = second->GetDynamicParticle();
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
             if(aSec->GetDefinition()==G4Neutron::Neutron()) nNeutrons++;
             if(isec==0) 
             {
               G4cout << aFinalState->GetNumberOfSecondaries();
             }
             else
             {
               G4cout << -1;
             }
             G4cout <<endl;
             delete second;
           }
           cout << "total secondaries"<<aFinalState->GetNumberOfSecondaries()<<" "<<
                   nNeutrons<<" "<< endl;
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   return EXIT_SUCCESS;
}


