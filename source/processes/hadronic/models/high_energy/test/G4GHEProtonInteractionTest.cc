// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GHEProtonInteractionTest.cc,v 1.1 1999-01-08 16:33:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>

#include "G4Material.hh"

#include "G4ProcessManager.hh"
#include "G4GHEProtonInelastic.hh"
#include "G4GHEProtonInelasticProcess.hh"

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

int main()
{
   G4cout.setf( ios::scientific, ios::floatfield );
   ofstream outFile( "ProtonInelastic.listing", ios::out);
   outFile.setf( ios::scientific, ios::floatfield );

  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3 );
  G4Material* Fe = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3 );
  G4Material* Pb = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3 );

  G4Element*   H = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Material* wa = new G4Material ("Water" , 1.*g/cm3, 2);
  wa->AddElement(H,2);
  wa->AddElement(O,1);

//  G4Element::DumpInfo(); 
//  G4Material::DumpInfo();

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int imat = 0;

    G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);

    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                  (*theMaterialTable)(imat),
						  "LFrame", 0, 0, 0);

    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                          "PFrame",LogicalFrame,0,false,0);

  G4ParticleDefinition* Proton = G4Proton::ProtonDefinition();

  //--------- Processes definition ---------
  G4ProcessManager* theProtonProcessManager = Proton->GetProcessManager();

  G4GHEProtonInelasticProcess theProtonInelasticProcess; 
  theProtonInelasticProcess.SetDEBUG( 1 );
  G4GHEProtonInelastic theProtonInelastic;
  theProtonInelastic.SetDEBUG( 1 );
  theProtonInelasticProcess.RegisterMe(&theProtonInelastic);
  theProtonProcessManager->AddProcess(&theProtonInelasticProcess);
  G4ForceCondition* condition = NULL;

  // ------- Build CrossSection Tables -------
  Proton->SetCuts(1.*mm);

  G4double protonEnergy = 100.*GeV;
  G4ParticleMomentum protonDirection(0.,0.,1.);
  G4DynamicParticle aProton(G4Proton::Proton(),protonDirection,protonEnergy);
  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;

  G4Track* ptrack = new G4Track(&aProton,aTime,aPosition) ;
  G4Track& aTrack = (*ptrack) ;
  ptrack->SetVolume(PhysicalFrame);

  G4Step* Step = new G4Step();  G4Step& aStep = (*Step);
  Step->SetTrack(ptrack);
 
  G4int nkin = 10; 
  // --------- Test the GetMeanFreePath                   
  G4double meanFreePath;
  for ( G4int i=0 ; i<nkin ; i++)
    {
       aProton.SetKineticEnergy(100.*GeV);

       meanFreePath = theProtonInelasticProcess.GetMeanFreePath(aTrack, 1., condition);

       outFile <<"  " << "100 GeV"  << "            " << meanFreePath/cm << endl ;
    }
  outFile << " " << endl;
   
  // --------- Test the PostStepDoIt

    G4Material * aNewMaterial = (*theMaterialTable)(imat) ;
    LogicalFrame->SetMaterial(aNewMaterial); 
    aProton.SetKineticEnergy(100.*GeV);
    aProton.SetMomentumDirection(0., 0., 1.);
    G4VParticleChange* aFinalState; 
    G4double totalCrossSection   = 0.;
    G4double elasticCrossSection = 0.;
    G4Material* targetMaterial = aTrack.GetMaterial();
    const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
    G4Element* targetElement = theProtonInelasticProcess.GetElement(incidentParticle, targetMaterial,
                                                                    totalCrossSection, elasticCrossSection);
    printf("ProtonInelasticProcess: CrossSections %10.1f %10.1f \n",
                                    totalCrossSection/millibarn, elasticCrossSection/millibarn); 
    G4int iteration = 0;   
    do {
      // The PostStep- Method does not work so far
      //        aFinalState = theProtonInelasticProcess.PostStepDoIt(aTrack, aStep);


      printf("ProtonInelasticProcess: call ProduceYourself \n");
                                            
      G4VParticleChange* aFinalState = new G4VParticleChange();
       
      aFinalState = theProtonInelastic.ProduceYourself( aTrack, targetElement,
                                                             totalCrossSection, elasticCrossSection );

      //      printf("ProtonInelasticProcess: ResetNumberOfInteractionLengthLeft \n");
 
      //      theProtonInelasticProcess.ResetNumberOfInteractionLengthLeft();

         outFile << " -----------------------------------------------------------" << endl;   
        // check the kinematic
        //
        G4double Tkin0 = aProton.GetKineticEnergy();
        G4double Px0   = aProton.GetMomentum().x() ,
                 Py0   = aProton.GetMomentum().y() ,
                 Pz0   = aProton.GetMomentum().z() ;
        outFile << " Initial Proton:  Px= " << Px0 << "   Py= " << Py0 << "   Pz= " << Pz0 
                           << "   Tkin= " << Tkin0 << endl ;

        G4double Tkin1 =  aFinalState->GetEnergyChange();
        G4double Px1   = (aFinalState->GetMomentumChange()->x())*Tkin1 ,

                 Py1   = (aFinalState->GetMomentumChange()->y())*Tkin1 ,
                 Pz1   = (aFinalState->GetMomentumChange()->z())*Tkin1 ;
        G4double Edep  =  aFinalState->GetLocalEnergyDeposited();
        outFile << "   final State:  Px= " << Px1 << "   Py= " << Py1 << "   Pz= " << Pz1 
                           << "   Tkin= " << Tkin1 << "   Edep= " << Edep << endl ;

        G4double Tkin2 = 0., Px2 = 0., Py2 = 0., Pz2 = 0. ;
        G4Track * aSecondary;
        for (G4int i=0; i<aFinalState->GetNumberOfSecondaries(); i++)
            {
            aSecondary = aFinalState->GetSecondary(i) ;
            Tkin2 =  aSecondary->GetKineticEnergy() ;
            Px2   = (aSecondary->GetMomentum()).x() ;
            Py2   = (aSecondary->GetMomentum()).y() ;
            Pz2   = (aSecondary->GetMomentum()).z() ;
            outFile << "  secondary: Px= " << Px2 << "   Py= " << Py2 << "   Pz= " << Pz2
                           << "   Tkin= " << Tkin2 << endl ;
           // NOTE - Secondaries are normally deleted by the track destructor !
	    //            delete aFinalState->GetSecondary(i);
           }
	//        aFinalState->Clear();

        iteration++; 

   }  while (iteration < 10) ;

  return EXIT_SUCCESS;
}
