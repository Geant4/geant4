// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GHETest.cc,v 1.2 1999-12-15 14:52:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"

#include "G4ProcessManager.hh"
#include "G4GHEProtonInelastic.hh"
#include "G4GHEProtonInelasticProcess.hh"

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"

int main()
{

  G4ParticleDefinition* Proton = G4Proton::ProtonDefinition();

  //--------- Processes definition ---------
  G4ProcessManager* theProtonProcessManager = Proton->GetProcessManager();

  G4GHEProtonInelasticProcess theProtonInelasticProcess; 
  theProtonInelasticProcess.SetDEBUG( 1 );
  G4GHEProtonInelastic theProtonInelastic;
  theProtonInelastic.SetDEBUG( 1 );
  theProtonInelastic.SetMaxNumberOfSecondaries(128);
  theProtonInelasticProcess.RegisterMe(&theProtonInelastic);
  theProtonProcessManager->AddProcess(&theProtonInelasticProcess);
  G4ForceCondition* condition = NULL;


  G4Element* targetElement = new G4Element("Iron", "Fe", 26., 55.85*g/mole);
  G4double protonEnergy = 100.*GeV;
  G4ParticleMomentum protonDirection(0.,0.,1.);
  G4DynamicParticle aProton(G4PionPlus::PionPlus(),protonDirection,protonEnergy);
  G4DynamicParticle aPion(G4Proton::Proton(),protonDirection,protonEnergy);
  aProton.SetDefinition( aPion.GetDefinition() ); 
  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0.;
  G4Track* ptrack = new G4Track(&aProton,aTime,aPosition);
  G4Track& aTrack = (*ptrack); 


  G4VParticleChange* aFinalState; 
  G4double totalCrossSection   = 1.;
  G4double elasticCrossSection = 0.;
  G4int iteration = 0;   
    do {

      printf("ProtonInelasticProcess: call ProduceYourself \n");
                                            
      G4VParticleChange* aFinalState = new G4VParticleChange();
       
      aFinalState = theProtonInelastic.ProduceYourself( aTrack, targetElement,
                                                             totalCrossSection, elasticCrossSection );

      G4Track* aSecondary;
      G4double Tkin, Px,Py,Pz;

      for (G4int i=0; i<aFinalState->GetNumberOfSecondaries(); i++)
          {
            aSecondary = aFinalState->GetSecondary(i) ;
            Tkin  =  aSecondary->GetKineticEnergy() ;
            Px    = (aSecondary->GetMomentum()).x() ;
            Py    = (aSecondary->GetMomentum()).y() ;
            Pz    = (aSecondary->GetMomentum()).z() ;
            printf("aSecondary: %10.1f %10.1f %10.1f %10.1f \n",
                    Tkin,Px,Py,Pz);  


           // NOTE - Secondaries are normally deleted by the track destructor !
	    //            delete aFinalState->GetSecondary(i);
           }
	//        aFinalState->Clear();

        iteration++; 

   }  while (iteration < 10) ;

  return EXIT_SUCCESS;
}
