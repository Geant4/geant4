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
// $Id: G4IonTest.cc,v 1.1 2003-10-08 12:19:47 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
 
#include "G4IonInelasticProcess.hh"
#include "G4LowEIonFragmentation.hh"

#include "G4DynamicParticle.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

int main()
  {
    G4cout.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile( "InInelasticAlpha.listing.GetMeanFreePath", std::ios::out);
    outFile.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile1( "InInelasticAlpha.listing.DoIt", std::ios::out);
    outFile1.setf( std::ios::scientific, std::ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
 // constructing the particles
 
 G4LeptonConstructor aC1;
 G4BaryonConstructor aC2;
 G4MesonConstructor aC3;
 G4IonConstructor aC4;
 
 aC1.ConstructParticle();
 aC2.ConstructParticle();
 aC3.ConstructParticle();
 aC4.ConstructParticle();
    
    // ----------- now get all particles of interest ---------
   G4int numberOfParticles = 1;
   G4ParticleDefinition* theParticles[1];
   G4ParticleDefinition* theIon = G4ParticleTable::GetParticleTable()->FindIon(6, 12, 0, 6);
   theParticles[0]=theIon;
   
      
   G4IonInelasticProcess theInelasticProcess; 
   G4LowEIonFragmentation theModel;
   G4cout << "Model instanciated!!!"<<G4endl;
   
   G4cout << "Done with all the process definitions"<<G4endl;
   //   G4cout << "Building the CrossSection Tables. This will take a while"<<G4endl;
   
   // define the Target, U in this case.
   G4Nucleus theTarget(16, 8);
   
   // ----------- needed to Build a track -------------------
   
   G4ParticleMomentum theDirection(0.,0.,1.);
   G4ThreeVector aPosition(0.,0.,0.);
   G4double aTime = 0.0;
   
   G4double incomingEnergy;
   G4int k, i, hpw = 0;   
   // --------- Test the ApplyYourself --------------
   G4cout << "Entering the Apply loops!!!!!"<< G4endl;
   G4ParticleChange* aFinalState;
   G4cout << "Test the Apply: please enter the number of events"<<G4endl;
   G4int ll0;
   G4cin >> ll0;
   G4int debugThisOne=1;
   G4cout << "Please enter the neutron energy"<<G4endl;
   G4cin >> incomingEnergy;
         for( G4int l=0; l<ll0; ++l )
         {
           G4cerr << "Event number "<<l<<endl;
	   G4DynamicParticle* aParticle =
             new G4DynamicParticle( theParticles[i], theDirection, incomingEnergy );
           G4Track* aTrack = new G4Track( aParticle, aTime, aPosition );
           ++hpw;
           if(hpw == 1000*(hpw/1000))
           G4cerr << "FINAL EVENTCOUNTER=" << hpw
                << " current energy: " << incomingEnergy
                << " of particle " << aParticle->GetDefinition()->GetParticleName() << G4endl;
           if (hpw==debugThisOne)
           {
            debugThisOne+=0;
           }
           G4cout << "Last chance before DoIt call: "
                <<G4endl;
           aFinalState = (G4ParticleChange*)  (theModel.ApplyYourself( *aTrack,  theTarget));
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
             second = aFinalState->GetSecondary(isec);
             aSec = second->GetDynamicParticle();
             G4cout << "SECONDARIES info";
             G4cout << aSec->GetTotalEnergy();
             G4cout << aSec->GetMomentum();
	     G4cout << (1-isec)*aFinalState->GetNumberOfSecondaries();
	     G4cout << G4endl;
	     G4ParticleDefinition * currentDefinition = second->GetDefinition();
	     G4cout << "A and Z of secondary "<<currentDefinition->GetBaryonNumber()<<" "<<currentDefinition->GetPDGCharge()<<G4endl;
             delete second;
           }
           delete aParticle;
           delete aTrack;
           aFinalState->Clear();
         }  // event loop
   return EXIT_SUCCESS;
}


