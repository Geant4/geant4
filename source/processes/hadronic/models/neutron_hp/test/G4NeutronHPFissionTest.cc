// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFissionTest.cc,v 1.2 1999-12-15 14:53:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "../src/G4NeutronHPIsoData.cc"
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4HadronFissionProcess.hh"
 
#include "G4LFission.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPIsoData.hh"

#include "G4DynamicParticle.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

#include "g4templates.hh"
#include "G4NeutronHPChannel.hh"
#include "G4NeutronHPElasticFS.hh"
#include "G4NeutronHPFissionFS.hh"
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

//    Natural Carbon moved to: 6_12_Carbon
     G4Material *thePS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
     G4Element *elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
     G4Element *elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
     thePS->AddElement( elC, 8 );
     thePS->AddElement( elH, 8 );
     theMaterials[1] = thePS;

// dead loop or out of memory between event 2000 and 3000
     G4Material *theLi = new G4Material(name="Lithium", density=0.534*g/cm3, nEl=1);
     G4Element *elLi = new G4Element(name="Lithium", symbol="Li", iz=3., a=6.941*g/mole);
     theLi->AddElement( elLi, 1);
     theMaterials[2] =theLi ;
     

// dead loop or out of memory between event 5000 and 6000
     G4Material *theB = new G4Material(name="Boron", density=0.9*g/cm3, nEl=1);
     G4Element *elB = new G4Element(name="Boron", symbol="B", iz=5, a=10.811*g/mole);
     theB->AddElement( elB, 1);
     theMaterials[3] = theB;

     G4Material *theN = new G4Material(name="Nitrogen", density=0.9*g/cm3, nEl=1);
     G4Element *elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
     theN->AddElement( elN, 1 );
     theMaterials[4] = theN;

    // redundant energy distributions not
    // given; FS/8_16_Oxygen reaches EOF at G4NeutronHPFissionFS.cc line 73
     G4Material *theO = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
     G4Element *elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
     theO->AddElement( elO,  1 );
     theMaterials[5] = theO;

//     // no FS data     
//     G4Material *theF = new G4Material(name="Fluor", density=1.1*g/cm3, nEl=1);
//     G4Element *elF = new G4Element(name="Fluor", symbol="F", iz=9, a=18.9984*g/mole);
//     theF->AddElement( elF,  1 );
//     theMaterials[6] = theF;
 
// Natural moved to 12_24_Magnesium
     G4Material *theMg = new G4Material(name="Magnesium", density=10.*g/cm3, nEl=1);
     G4Element *elMg = new G4Element(name="Magnesium", symbol="Mg", iz=12, a=24.305*g/mole);
     theMg->AddElement( elMg,  1 );
     theMaterials[7] = theMg;
     
     G4Material *theAl = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
     G4Element *elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
     theAl->AddElement( elAl, 1 );
     theMaterials[8] = theAl;
     
// Natural moved to 14_28_Silicon
    // redundant energy distributions not
    // given; FS/14_28_Silicon reaches EOF at G4NeutronHPFissionFS.cc line 73
     G4Material *theSi = new G4Material(name="Silicon", density=2.33*g/cm3, nEl=1);
     G4Element *elSi  = new G4Element(name="Silicon", symbol="Al", iz=14., a=28.0855*g/mole);
     theSi->AddElement( elSi, 1 );
     theMaterials[9] = theSi;
     
// Natural moved to 17_35_Chlorine
    // redundant energy distributions not
    // given; FS/17_35_Chlorine reaches EOF at G4NeutronHPFissionFS.cc line 73
     G4Material *theCl = new G4Material(name="Chlorine", density=3*g/cm3, nEl=1);
     G4Element *elCl  = new G4Element(name="Chlorine", symbol="Al", iz=17, a=35.4527*g/mole);
     theCl->AddElement( elCl, 1 );
     theMaterials[10] = theCl;
     
// // only cross-sections, no FS data in END/B-VI
//      G4Material *theLAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
//      G4Element *elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
//      theLAr->AddElement( elAr, 1 );
//     theMaterials[11] = theLAr;
    
// Natural moved to 20_40_Calcium
     G4Material *theCa= new G4Material(name="Calcium", density=5*g/cm3, nEl=1);
     G4Element *elCa  = new G4Element(name="Calcium", symbol="Ca", iz=20., a=40.078*g/mole);
     theCa->AddElement( elCa, 1 );
     theMaterials[12] = theCa;
    
     G4Material *theCr= new G4Material(name="Chromium", density=1.393*g/cm3, nEl=1);
     G4Element *elCr  = new G4Element(name="Chromium", symbol="Cr", iz=24., a=39.95*g/mole);
     theCr->AddElement( elCr, 1 );
     theMaterials[13] = theCr;

     G4Material *theFe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
     G4Element *elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
     theFe->AddElement( elFe, 1 );
     theMaterials[14] = theFe;

     G4Material *theNi = new G4Material(name="Nickel", density=7.87*g/cm3, nEl=1);
     G4Element *elNi = new G4Element(name="Nickel", symbol="Ni", iz=28., a=58.6934*g/mole);
     theNi->AddElement( elNi, 1 );
     theMaterials[15] = theNi;

     G4Material *theCo = new G4Material(name="Cobalt", density=7.87*g/cm3, nEl=1);
     G4Element *elCo = new G4Element(name="Cobalt", symbol="Co", iz=27., a=58.9332*g/mole);
     theCo->AddElement( elCo, 1 );
     theMaterials[16] =theCo ;

     G4Material *theMb = new G4Material(name="Molybden", density=7.87*g/cm3, nEl=1);
     G4Element *elMb = new G4Element(name="Molybden", symbol="Mo", iz=42., a=95.94*g/mole);
     theMb->AddElement( elMb, 1 );
     theMaterials[17] = theMb;

     G4Material *theBa = new G4Material(name="Barium", density=7.87*g/cm3, nEl=1);
     G4Element *elBa = new G4Element(name="Barium", symbol="Ba", iz=56., a=137.327*g/mole);
     theBa->AddElement( elBa, 1 );
     theMaterials[18] = theBa;

     G4Material *theTa = new G4Material(name="Tantalum", density=7.87*g/cm3, nEl=1);
     G4Element *elTa = new G4Element(name="Tantalum", symbol="Ta", iz=73., a=180.9479*g/mole);
     theTa->AddElement( elTa, 1 );
     theMaterials[19] = theTa;

     G4Material *theW  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
     G4Element *elW  = new G4Element(name="Tungston", symbol="W", iz=74., a=183.85*g/mole);
     theW->AddElement( elW, 1 );
     theMaterials[20] = theW;

      G4Material *thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
      G4Element *elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
      thePb->AddElement( elPb, 1 );
     theMaterials[21] = thePb;

     G4Material *theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
      G4Element *elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
      theU->AddElement( elU, 1 );
     theMaterials[22] = theU;
    
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
   G4HadronFissionProcess* theProcesses[1];
      
   G4ProcessManager* theNeutronProcessManager = new G4ProcessManager(theNeutron);
   theNeutron->SetProcessManager(theNeutronProcessManager);
   G4HadronFissionProcess theFissionProcess; 
   G4LFission theFission;
   G4NeutronHPFission theNeutronHPFission;
   G4cout << "Fission instanciated!!!"<<G4endl;
   theFissionProcess.RegisterMe(&theFission);
   theFission.SetMinEnergy(21*MeV);
   theFissionProcess.RegisterMe(&theNeutronHPFission);
   theNeutronProcessManager->AddDiscreteProcess(&theFissionProcess);
   theProcesses[0] = &theFissionProcess;
   
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
   G4cout <<"Now debug the DoIt: enter the problem event number"<< G4endl;
   G4int debugThisOne=1;
   G4cin >> debugThisOne;
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
           if(hpw==1000*(hpw/1000))
           G4cerr << "FINAL EVENTCOUNTER=" << hpw
                << " current energy: " << incomingEnergy
                << " of particle " << aParticle->GetDefinition()->GetParticleName() 
                << " in material " << theMaterials[k]->GetName() << G4endl;
           if (hpw==debugThisOne)
           {
             debugThisOne+=0;
           }
           G4cout << "Last chance before DoIt call: "
//                << theNeutronHPFission.GetNiso()
                <<G4endl;
           aFinalState =
	   (G4ParticleChange*) (theProcesses[i]->PostStepDoIt( *aTrack, aStep ));
           G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries();
           G4double theFSEnergy = aFinalState->GetEnergyChange();
           G4ThreeVector * theFSMomentum= aFinalState->GetMomentumChange();
           G4cout << "FINAL STATE = "<<theFSEnergy<<" ";
           G4cout <<*theFSMomentum<<G4endl;
           G4Track * second;
           G4DynamicParticle * aSec;
           G4int isec;
           for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
           {
             G4cout << "SECONDARIES info";
             second = aFinalState->GetSecondary(isec);
             aSec = second->GetDynamicParticle();
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
             if(isec==0) 
             {
               G4cout << aFinalState->GetNumberOfSecondaries();
             }
             else
             {
               G4cout << -1;
             }
             G4cout <<G4endl;
             delete second;
           }
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   return EXIT_SUCCESS;
}


