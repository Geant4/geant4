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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Timer.hh"
 
#include "G4Material.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
 
#include "G4ProcessManager.hh"
#include "G4HadronFissionProcess.hh"
 
#include "G4ParaFissionModel.hh"

#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4Step.hh"

#include "NametoGeant3Number.cc"
 
 int main()
  {
    //    G4cout.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile( "InelasticAlpha.listing.GetMeanFreePath", std::ios::out);
    outFile.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile1( "InelasticAlpha.listing.DoIt", std::ios::out);
    outFile1.setf( std::ios::scientific, std::ios::floatfield );
    
    RanecuEngine theEngine;
    HepRandom::setTheEngine( &theEngine );
    theEngine.showStatus();
    G4cout << G4endl;
    
    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
    G4Material *theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
    G4Element *elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
    theU->AddElement( elU, 1 );
    
    G4int numberOfMaterials = 1;
    G4Material *theMaterials[5];
    theMaterials[0] = theU;
            
    // ----------- the following is needed for building a track ------------
    
    static const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
    G4int imat = 0;
    G4Box* theFrame = new G4Box ( "Frame",10*m, 10*m, 10*m );
    
    G4LogicalVolume *LogicalFrame = new G4LogicalVolume( theFrame,
                                                         theMaterials[0],
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

    G4BaryonConstructor Barions;
    Barions.ConstructParticle();

        G4ShortLivedConstructor ShortLived;
        ShortLived.ConstructParticle();
 
    G4int numberOfParticles = 1;
    G4ParticleDefinition *theParticles[5];
    G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
    theParticles[ 0] = theNeutron;
    
    //------ all the particles are Done ----------
    //------ Processes definitions Follow ---------

    // this will be the model class for fission
    G4ParaFissionModel * theModel = new G4ParaFissionModel;
    
    G4HadronFissionProcess *theProcesses[1];
    
    G4ProcessManager *theNeutronProcessManager = new G4ProcessManager(theNeutron);
    theNeutron->SetProcessManager(theNeutronProcessManager);
    G4HadronFissionProcess *theNeutronFissionProcess =
      new G4HadronFissionProcess(); 
    theNeutronFissionProcess->RegisterMe( theModel );
    theNeutronProcessManager->AddDiscreteProcess( theNeutronFissionProcess );
    theProcesses[0] = theNeutronFissionProcess;
        
    G4ForceCondition *condition = new G4ForceCondition;
    *condition = NotForced;
    
    G4ParticleMomentum theDirection( 0., 0., 1. );
    G4ThreeVector aPosition( 0., 0., 0. );
    G4double aTime = 0.0;
    
    G4Step aStep;
    G4double meanFreePath;
    G4double incomingEnergy;
    G4ParticleDefinition *currentDefinition;
    G4Track *currentTrack;
    
    G4int kl, kr;
      kl = 0;
      kr = kl;
    
    G4int il, ir;
      il = 0;
      ir = il;


    G4cout << "Please enter the initial kinetic energy (MeV): "<< std::flush;
    G4cin >> incomingEnergy;
    incomingEnergy *= MeV;
    
    G4int nEvents;
    G4cout << "Please enter the number of events: "<< std::flush;
    G4cin >> nEvents;
    
    FILE *g4data = NULL;
    if( kl == kr && il == ir ) {
      if( (g4data = fopen("g4data","w")) == NULL )
	G4cerr << "Cannot open g4data file" << G4endl;
    }
    if( !g4data )
      G4cout << "Output to file will be forbidden!"<< G4endl;
   
    G4int nCPU = 0;
    G4double dCPU = 0.0;
    G4Timer timerEvent;
    G4Timer timerTotal;
    timerTotal.Start();
 
    G4int k, i;
    for( k = kl; k <= kr; k++ )
      for( i = il; i <= ir; i++ ) {
	theParticles[i]->SetCuts( (i ? 1.0 : 0.01)*mm );
	LogicalFrame->SetMaterial( theMaterials[0] ); 
    
    // --------- Test the PostStepDoIt now  --------------
    
	G4ParticleChange *aFinalState;
//	G4cout << "debugging the material "<<k<<" "<<theMaterials[0]<<" end test"<<G4endl;
        G4int it = theMaterials[0]; G4cout << it<<endl;
	
	G4DynamicParticle aParticle( theParticles[i], theDirection, incomingEnergy );
	if( g4data ) {
	  fprintf( g4data, "%s\n", theParticles[i]->GetParticleName() );
	  fprintf( g4data, "%s\n", theMaterials[0]->GetName() );
	  fprintf( g4data, "%6.1f\n", incomingEnergy/GeV, nEvents );
	}
	G4double currentPx, currentPy, currentPz;
	G4double kineticEnergy, totalEnergy, currentMass;
    
	G4int eventCounter = 0;
	G4int l;
	for( l=0; l<nEvents; ++l ) {
	  aParticle.SetDefinition( theParticles[i] );
	  aParticle.SetMomentumDirection( theDirection );
	  aParticle.SetKineticEnergy( incomingEnergy );
      
	  G4Track *aTrack = new G4Track( &aParticle, aTime, aPosition );
	  aTrack->SetTouchable( theTouchable );
	  aTrack->SetStep( &aStep );
	  aStep.SetTrack( aTrack );
	  G4cout << "Event:" << std::setw(4) << l+1 << 
	            " Material:" << theMaterials[0]->GetName() <<
	            " Particle:" << theParticles[i]->GetParticleName() << '\t' << std::flush;

	  timerEvent.Start();
	  aFinalState = (G4ParticleChange * ) (theProcesses[i]->PostStepDoIt( *aTrack, aStep ) );
	  timerEvent.Stop();

	  G4int nSec = aFinalState->GetNumberOfSecondaries();
	  //delete aTrack;
      
	  // prepare the analysis current quantities
      
	  G4cout << "Definitions:";
	  G4int istart, ipar;        
	  if( aFinalState->GetStatusChange() == fAlive ) {// current particle is still alive
	    istart = 0;
	    G4cout << '1';
	    currentDefinition = aParticle.GetDefinition();
	    kineticEnergy = aFinalState->GetEnergyChange()/GeV;
	    currentMass = currentDefinition->GetPDGMass()/GeV;
	    totalEnergy = kineticEnergy+currentMass;
	    const G4ParticleMomentum *mom = aFinalState->GetMomentumChange();
	    G4double p = sqrt( kineticEnergy*kineticEnergy + 2.0*kineticEnergy*currentMass );
	    currentPx = mom->x()*p;
	    currentPy = mom->y()*p;
	    currentPz = mom->z()*p;
	    ipar = NameToGeant3Number( currentDefinition );
	    if( g4data )
	      fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d %13.6E %s\n",
		       kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
		       ipar, nSec+1, l, currentMass, currentDefinition->GetParticleName() );
	  } else if (nSec !=0) {
	    istart = 1;
	    currentTrack = aFinalState->GetSecondary(0);
	    G4cout << '2';
	    currentDefinition = currentTrack->GetDefinition();
	    G4cout << "A and Z of secondary "<<currentDefinition->GetBaryonNumber()<<" "<<currentDefinition->GetPDGCharge()<<G4endl;
	    kineticEnergy = currentTrack->GetKineticEnergy()/GeV;
	    currentMass = currentDefinition->GetPDGMass()/GeV;
	    totalEnergy = kineticEnergy+currentMass;
	    currentPx = currentTrack->GetMomentum().x()/GeV;
	    currentPy = currentTrack->GetMomentum().y()/GeV;
	    currentPz = currentTrack->GetMomentum().z()/GeV;
	    ipar = NameToGeant3Number( currentDefinition );
	    if( g4data )
	      fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d %13.6E %s\n",
		       kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
		       ipar, nSec, l, currentMass, currentDefinition->GetParticleName() );
	  }
	  G4cout << "/3:" << nSec - istart << 
	    " CPU Time(s):" << std::setw(5) << (dCPU += timerEvent.GetSystemElapsed()) << '/' << ++nCPU << G4endl;
	  for( G4int ii=istart; ii<nSec; ++ii ) {
	    G4Track *currentTrack; 
	    currentTrack = aFinalState->GetSecondary(ii);
	    currentDefinition = currentTrack->GetDefinition();
	    kineticEnergy = currentTrack->GetKineticEnergy()/GeV;
	    currentMass = currentDefinition->GetPDGMass()/GeV;
	    totalEnergy = kineticEnergy+currentMass;
	    currentPx = currentTrack->GetMomentum().x();
	    currentPy = currentTrack->GetMomentum().y();
	    currentPz = currentTrack->GetMomentum().z();
	    ipar = NameToGeant3Number( currentDefinition );
	    if( g4data )
	      fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d %13.6E %s\n",
		       kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
		       ipar, 0, l, currentMass, currentDefinition->GetParticleName() );
	  }
	  aFinalState->Clear();
	} // end of event loop
      } // end of i loop
    timerTotal.Stop();    
    G4cout << "\tTotal time " << timerTotal <<
              "\n\tPure function time:" << dCPU <<
              "\n\tTotal number of events:" << nCPU << G4endl; 
    fclose( g4data );
    return EXIT_SUCCESS;
}
 /* end of code */
 
