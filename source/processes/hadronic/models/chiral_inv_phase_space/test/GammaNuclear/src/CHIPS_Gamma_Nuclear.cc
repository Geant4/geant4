//#define debug
#define edebug
    
#include "GammaNuclear/src/ParticleInfo.h"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Timer.hh"
 
#include "G4Material.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
 
#include "G4ProcessManager.hh"
 
#include "G4GammaNuclearReaction.hh"

#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"

int main()
{
    
  RanecuEngine theEngine;
  HepRandom::setTheEngine( &theEngine );
  theEngine.showStatus();
  G4cout << G4endl;
  
  G4String name, symbol;
  G4double a, iz, density;
  G4int nEl;
  
  G4Material *theCu = new G4Material(name="Copper", density=8.96*g/cm3, nEl=1);
  G4Element *elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a=63.55*g/mole);
  theCu->AddElement( elCu, 1 );
  G4Material *thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
  G4Element *elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
  thePb->AddElement( elPb, 1 );
  G4Material *theFe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
  G4Element *elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
  theFe->AddElement( elFe, 1 );
  G4Material *theW  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
  G4Element *elW  = new G4Element(name="Tungston", symbol="W", iz=74., a=183.85*g/mole);
  theW->AddElement( elW, 1 );
  G4Material *theLAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
  G4Element *elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
  theLAr->AddElement( elAr, 1 );
  G4Material *thePS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
  G4Element *elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
  G4Element *elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
  thePS->AddElement( elC, 8 );
  thePS->AddElement( elH, 8 );
  G4Material *thePbWO4 = new G4Material(name="PbWO4", density=12.0*g/cm3, nEl=3);
  // approximative number
  G4Element *elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
  thePbWO4->AddElement( elPb, 1 );
  thePbWO4->AddElement( elW,  1 );
  thePbWO4->AddElement( elO,  4 );
  // approximate numbers for O
  G4Material *theO = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
  theO->AddElement( elO,  1 );
  G4Material *theBe = new G4Material(name="Beryllium", density=1.848*g/cm3, nEl=1);
  G4Element *elBe  = new G4Element(name="Beryllium", symbol="Be", iz=4., a=9.01*g/mole);
  theBe->AddElement( elBe, 1 );
  G4Material *theAl = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
  G4Element *elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
  theAl->AddElement( elAl, 1 );
  G4Material *theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
  G4Element *elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
  theU->AddElement( elU, 1 );
  G4Material *theBGO = new G4Material(name="BGO", density=2.15*g/cm3, nEl=3);
  G4Element *elBi = new G4Element(name="Bismuth", symbol="Bi", iz=83., a=208.98*g/mole);
  G4Element *elGe = new G4Element(name="Germanium", symbol="Ge", iz=32., a=72.59*g/mole);
  theBGO->AddElement( elBi, 4 );
  theBGO->AddElement( elGe, 3 );
  theBGO->AddElement( elO, 12 );
  G4Material *theNaI = new G4Material(name="NaI", density=3.67*g/cm3, nEl=2);
  G4Element *elNa = new G4Element(name="Sodium", symbol="Na", iz=11., a=22.990*g/mole);
  G4Element *elI = new G4Element(name="Iodine", symbol="I", iz=53., a=126.904*g/mole);
  theNaI->AddElement( elNa, 1 );
  theNaI->AddElement( elI, 1 );
  G4Material *theCsI = new G4Material(name="CsI", density=4.53*g/cm3, nEl=2);
  G4Element *elCs = new G4Element(name="Cesium", symbol="Cs", iz=55., a=132.905*g/mole);
  theCsI->AddElement( elCs, 1 );
  theCsI->AddElement( elI, 1 );
  G4Material *theKapton = new G4Material(name="Kapton", density=1.53*g/cm3, nEl=4); 
  // formula: private communications, see mail.
  theKapton->AddElement( elC, 22 );
  theKapton->AddElement( elH, 10 );
  theKapton->AddElement( elO,  5 );
  G4Element *elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
  theKapton->AddElement( elN, 2 );
  G4Material *theH = new G4Material(name="Hydrogen", density=1.53*g/cm3, nEl=1); 
  theH->AddElement( elH, 1 );
  G4Material *theC = new G4Material(name="Carbon", density=1.032*g/cm3, nEl=1);
  theC->AddElement( elC, 1 );
  G4Material *theLi6 = new G4Material(name="Li6", density=2.70*g/cm3, nEl=1);
  G4int nIso;
  G4Element *elLi6  = new G4Element(name="Li6", symbol="Li", nIso = 1);
  G4Isotope * isoLi6 = new G4Isotope(name="Li6", 3, 6, a=6.*g/mole);
  elLi6->AddIsotope(isoLi6, 1);
  theLi6->AddElement(elLi6 , 1 );
  G4Material *theTa = new G4Material(name="Tantalum", density=10.*g/cm3, nEl=1); // dens ? 
  G4Element *elTa = new G4Element(name="Tantalum", symbol="Ta", iz=73., a=181*g/mole);
  theTa->AddElement( elTa, 1 );
  G4Material *theCa = new G4Material(name="Calcium", density=6.*g/cm3, nEl=1); // dens ?
  G4Element *elCa = new G4Element(name="Calcium", symbol="Ca", iz=20., a=40*g/mole);
  theCa->AddElement( elCa, 1 );
  
  const G4int numberOfMaterials = 20;
  G4Material *theMaterials[numberOfMaterials];
  theMaterials[ 0] = theCu;
  theMaterials[ 1] = thePb;
  theMaterials[ 2] = theFe;
  theMaterials[ 3] = theW;
  theMaterials[ 4] = theLAr;
  theMaterials[ 5] = thePS;
  theMaterials[ 6] = thePbWO4;
  theMaterials[ 7] = theO;
  theMaterials[ 8] = theBe;
  theMaterials[ 9] = theAl;
  theMaterials[10] = theU;
  theMaterials[11] = theBGO;
  theMaterials[12] = theNaI;
  theMaterials[13] = theCsI;
  theMaterials[14] = theKapton;
  theMaterials[15] = theH;
  theMaterials[16] = theC;
  theMaterials[17] = theLi6;
  theMaterials[18] = theTa;
  theMaterials[19] = theCa;
 
  // ----------- here all materials have been defined -----------    
  // ----------- the following is needed for building a track ------------
  
  static const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
  G4int imat = 0;
  G4Box* theFrame = new G4Box ( "Frame",10*m, 10*m, 10*m );
  
  G4LogicalVolume *LogicalFrame = new G4LogicalVolume( theFrame,
                                                       (*theMaterialTable)[imat],
                                                       "LFrame", 0, 0, 0 );
  
  G4PVPlacement *PhysicalFrame = new G4PVPlacement( 0, G4ThreeVector(),
                                                   "PFrame", LogicalFrame, 0, false, 0 );
  
  G4RotationMatrix theNull;
  G4ThreeVector theCenter(0,0,0);
  G4GRSVolume * theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
  // ----------- now get all particles of interest ---------
 
  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BaryonConstructor Baryons;
  Baryons.ConstructParticle();

  G4IonConstructor Ions;
  Ions.ConstructParticle();

  G4ParticleDefinition *theGamma = G4Gamma::GammaDefinition();
  
  //------ all the particles are Done ----------
  //------ Model definitions Follow ---------

// @@@ add gamma nuclear reaction here.
  G4GammaNuclearReaction* theCHIPSGammaNuclear = new G4GammaNuclearReaction;
    
// some more stuff needed to create the track - to satisfy an internal consistency check

  G4ForceCondition *condition = new G4ForceCondition;
  *condition = NotForced;
  
  G4ParticleMomentum theDirection( 0., 0., 1. );
  G4ThreeVector aPosition( 0., 0., 0. );
  G4double aTime = 0.0;
  
  G4StepPoint aStepPoint;
  G4Step aStep;
  aStep.SetPreStepPoint(&aStepPoint);
  ///////////////////////////////////////////////G4double meanFreePath;
  G4double incomingEnergy;
  ///////////////////////////////////////////////G4ParticleDefinition *currentDefinition;
  ///////////////////////////////////////////////G4Track *currentTrack;

// end of stuff for the track
  G4int kl, kr;
  //do {
  //  G4cout << " 0) Copper      1) Lead         2) Iron      3) Tungsten" << G4endl;
  //  G4cout << " 4) LArgon      5) PolyStyrene  6) PbWO4     7) Oxygen" << G4endl;
  //  G4cout << " 8) Beryllium   9) Aluminium   10) Uranium  11) BGO" << G4endl;
  //  G4cout << "12) NaI        13) CsI         14) Kapton   15) Hydrogen" << G4endl;
  //  G4cout << "16) Carbon     17) 6_3_Lithium 18) Tantalum 19) Calcium" << G4endl;
  //  G4cout << "Please enter the material code" << G4endl;
  //  G4cout << "\tFrom: " << G4std::flush;
  //  G4cin >> kl;
  //  G4cout << "\tTo: " << G4std::flush;
  //  G4cin >> kr;
  //} while( kl < 0 || kr >= numberOfMaterials || kl > kr );
  // Now only for Carbon
  kl=16;
  kr=16;
  
  // energies of interest: 0.5, 1.0, 3.0, 12.0, 20.0 GeV

  G4cout << "Please enter the initial energy of gamma (in GeV): "<< G4std::flush;
  G4cin >> incomingEnergy;
  incomingEnergy *= GeV/MeV;
  
  G4int nEvents;
  G4cout << "Please enter the number of events: "<< G4std::flush;
  G4cin >> nEvents;
     
  G4int nCPU = 0;
  G4double dCPU = 0.0;
  G4Timer timerEvent;
  G4Timer timerTotal;
  timerTotal.Start();
    
// Prepare the analysis
  G4String file("../logs/liste.");
  G4String fileName; // p_al
  G4double weight; // sigma/A
  //weight = .058*millibarn; // Cross-section per nucleon @@ can be inside IF
  weight = 1.8*millibarn; // Cross-section @@ can be inside IF
  if(kl == 16)      fileName = "gamma_c";
  //else if(kl == 8)  fileName = "gamma_be";
  //else if(kl == 17) fileName = "gamma_li6";
  else if(kl == 9)  fileName = "gamma_al";
  //else if(kl == 0)  fileName = "gamma_cu";
  //else if(kl == 18) fileName = "gamma_ta";
  G4cout << "CHIPS_Gamma_Nuclear: filling "<<fileName<<" file"<<G4endl;
  G4String it;
  it = file+fileName;
  ANAParticleInfo theInformation(weight, it);
    
// Analysis prepared
 
  G4int k;
  for( k = kl; k <= kr; k++ )
  {
	LogicalFrame->SetMaterial( theMaterials[k] ); 
  
    // --------- Test the PostStepDoIt now  --------------
  
	G4ParticleChange *aFinalState;
	LogicalFrame->SetMaterial( theMaterials[k] ); 
	G4DynamicParticle aParticle( theGamma, theDirection, incomingEnergy );
	//G4double currentPx, currentPy, currentPz;
	//G4double kineticEnergy, totalEnergy, currentMass;
  
	//G4int eventCounter = 0;
	G4int l;
	for( l=0; l<nEvents; ++l ) {
	  aParticle.SetDefinition( theGamma );
	  aParticle.SetMomentumDirection( theDirection );
	  aParticle.SetKineticEnergy( incomingEnergy );
    
      G4Track *aTrack = new G4Track( &aParticle, aTime, aPosition );
	  aTrack->SetTouchable( theTouchable );
	  aTrack->SetStep( &aStep );
	  aStep.SetTrack( aTrack );
      aStepPoint.SetTouchable(theTouchable);
	  aStepPoint.SetMaterial(theMaterials[k]);
      aStep.SetPreStepPoint(&aStepPoint);
      aStep.SetPostStepPoint(&aStepPoint);
#ifdef debug
	  G4cout << "Event:" << G4std::setw(4) << l+1 << 
                " Material:" << theMaterials[k]->GetName() <<
	            " Particle:" << theGamma->GetParticleName() << '\t' << G4endl;
#endif
#ifdef edebug
	  if (!(l%100)) G4cout<<"Event #"<<l<<G4endl;
#endif
	  timerEvent.Start();
	  G4Nucleus theTarget(theMaterials[k]);
	  aFinalState=(G4ParticleChange*)(theCHIPSGammaNuclear->ApplyYourself(*aTrack,theTarget));
#ifdef debug
	  G4int nSec = aFinalState->GetNumberOfSecondaries();
      G4cout << "CHIPS_Gamma_Nuclear:  A NUMBER OF SECONDARIES = "<<nSec<<G4endl;
#endif
	  timerEvent.Stop();
	  // prepare the analysis current quantities
	  //G4int istart, ipar;  
	  //G4double etWithinCuts = 0.;    
	  //G4int ii;  
      G4Track* second;
	  G4DynamicParticle * aSec;
      G4int isec;
      for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
      {
        second = aFinalState->GetSecondary(isec);
        aSec = const_cast<G4DynamicParticle *> (second->GetDynamicParticle());
#ifdef debug
        G4cout << "SECONDARIES info Charge=";
        G4cout << aSec->GetDefinition()->GetPDGCharge()<<", (PDG or Z*1000+A)=";
	    if(aSec->GetDefinition()->GetPDGEncoding()) G4cout<<aSec->GetDefinition()->GetPDGEncoding();
	    else G4cout<<aSec->GetDefinition()->GetPDGCharge()*1000+aSec->
			   GetDefinition()->GetBaryonNumber();
        G4cout <<", E="<<aSec->GetTotalEnergy()<<", P="<<aSec->GetMomentum()<<G4endl;
#endif
	    // analysis part;
	    G4int theCHIPSCode = aSec->GetDefinition()->GetPDGEncoding();
	    if(theCHIPSCode == 0)
	      theCHIPSCode=static_cast<G4int>(aSec->GetDefinition()->GetPDGCharge()*1000+aSec->
                     GetDefinition()->GetBaryonNumber());
	    ANAParticle aPart(theCHIPSCode,
	                      aSec->GetMomentum().x(),
		                  aSec->GetMomentum().y(),
			              aSec->GetMomentum().z(),
			              aSec->GetTotalEnergy());
        theInformation.ProcessOne(aPart);
	    delete aSec;
	    delete second;
      }
      delete aTrack;
	  aFinalState->Clear();
	  if(l==500*(l/500)) 
	  {
	    theInformation.Analyse();
	    theInformation.Plot(fileName, l);
	  }
	} // End of event loop
  } // End of material loop
  timerTotal.Stop();    
  G4cout << "Terminating successfully"<<endl;
  G4cout << "\tTotal time " << timerTotal <<
            "\n\tPure function time:" << dCPU <<
            "\n\tTotal number of events:" << nCPU << G4endl; 
  return EXIT_SUCCESS;
}
 /* end of code */
 
