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
// $Id: G4LEppTest.cc,v 1.5 2001-10-11 08:54:02 fjones Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
#include "G4ProcessManager.hh"
#include "G4HadronElasticProcess.hh"
 
#include "G4LEpp.hh"

#include "G4DynamicParticle.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

#include "HbookHistogram.hh"


G4int sortEnergies(const double Px, const double Py, const double Pz,
                   const double Ekin, double* sortedPx, double* sortedPy,
                   double* sortedPz, double* sortedE)
{
   for(int i=0; i<10; ++i) {
      if(abs(Ekin) > sortedE[i]) {
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

   theHbookManager.SetFilename("G4LEppTest.hbook");

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
 
   G4int numberOfMaterials = 1;
   G4Material* theMaterials[2000];
   
   G4Element* elH = new G4Element(name="Hydrogen", symbol="H", 
                                  iz=1., a=1.01*g/mole);
   G4Material* theH = new G4Material(name="Hydrogen", density=1.032*g/cm3,
                                     nEl=1);
   theH->AddElement(elH, 1);
   theMaterials[1] = theH;
    
   G4int inputNumber = 1;
   theMaterials[0] = theMaterials[inputNumber];
   G4cout << "Active Material = " << theMaterials[0]->GetName() << G4endl;
    
   static const G4MaterialTable* theMaterialTable = 
      G4Material::GetMaterialTable();
   G4int imat = 0;   
   G4Box* theFrame = new G4Box("Frame", 10*m, 10*m, 10*m);
    
   G4LogicalVolume* LogicalFrame = 
      new G4LogicalVolume(theFrame, (*theMaterialTable)(imat),
                          "LFrame", 0, 0, 0);
    
   G4PVPlacement* PhysicalFrame = 
      new G4PVPlacement(0,G4ThreeVector(),
                        "PFrame", LogicalFrame, 0, false, 0);
   G4RotationMatrix theNull;
   G4ThreeVector theCenter(0, 0, 0);
   G4GRSVolume* theTouchable = 
      new G4GRSVolume(PhysicalFrame, &theNull, theCenter);

   // ----------- now get all particles of interest ---------
   G4int numberOfParticles = 1;
   G4ParticleDefinition* theParticles[1];
   G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
   theParticles[0] = theProton;
   
   //------ here all the particles are Done ----------
   G4cout << "Done with all the particles" << G4endl;
   G4cout << "Starting process definitions" << G4endl;

   //--------- Processes definitions ---------
   G4HadronElasticProcess* theProcesses[1];
      
   G4ProcessManager* theProtonProcessManager = new G4ProcessManager(theProton);
   theProton->SetProcessManager(theProtonProcessManager);
   G4HadronElasticProcess theProcess; 
   G4LEpp theModel;
   G4cout << "Model instanciated!!!" << G4endl;
   //   theModel.SetVerboseLevel(2);
   theProcess.RegisterMe(&theModel);
   theProtonProcessManager->AddDiscreteProcess(&theProcess);
   theProcesses[0] = &theProcess;
   
   G4ForceCondition* condition = new G4ForceCondition;
   *condition = NotForced;

   G4cout << "Done with all the process definitions" << G4endl;
   
   // ----------- needed to Build a Step, and a track -------------------
   
   G4ParticleMomentum theDirection(0., 1., 0.);
   //   G4double Pxyz = -sqrt(1./3.);
   //   G4ParticleMomentum theDirection(Pxyz, Pxyz, Pxyz);
   G4ThreeVector aPosition(0., 0., 0.);
   G4double aTime = 0.0;
   
   static G4Step aStep;
   static G4StepPoint aStepPoint;
   aStep.SetPreStepPoint(&aStepPoint);
   G4double meanFreePath;
   G4double incomingEnergy;

   G4ParticleChange* aParticleChange;

   G4cout << "Test the DoIt: please enter the number of events" << G4endl;
   G4int nEvent;
   G4cin >> nEvent;

   G4cout << "Please enter the proton energy in GeV" << G4endl;
   G4cin >> incomingEnergy;

   G4DynamicParticle* aParticle =
      new G4DynamicParticle(theParticles[0], theDirection, incomingEnergy*GeV);
   G4Track* aTrack = new G4Track(aParticle, aTime, aPosition);

   HbookHistogram* HPx1 = new HbookHistogram("Px1", 80, -1., 1.);
   HbookHistogram* HPy1 = new HbookHistogram("Py1", 80, -1., 1.);
   HbookHistogram* HPz1 = new HbookHistogram("Pz1", 80, -1., 1.);
   HbookHistogram* HPx2 = new HbookHistogram("Px2", 80, -1., 1.);
   HbookHistogram* HPy2 = new HbookHistogram("Py2", 80, -1., 1.);
   HbookHistogram* HPz2 = new HbookHistogram("Pz2", 80, -1., 1.);

   HbookHistogram* HAng1 = new HbookHistogram("Lab angle 1", 80, -180., 180.);
   HbookHistogram* HAng2 = new HbookHistogram("Lab angle 2", 80, -180., 180.);

   HbookHistogram* HNsec = 
           new HbookHistogram("Number of secondaries", 4, -0.5, 3.5);

   for (G4int i = 0; i < numberOfParticles; i++) {
      LogicalFrame->SetMaterial(theMaterials[0]); 
      for (G4int iEvent = 1; iEvent <= nEvent; iEvent++) {
         G4cout << "------------------------------ Event " << iEvent << G4endl;
         aTrack->SetTouchable(theTouchable);
         aStep.SetTrack(aTrack);
         aStepPoint.SetTouchable(theTouchable);
         aStepPoint.SetMaterial(theMaterials[0]);
         aStep.SetPreStepPoint(&aStepPoint);
         aStep.SetPostStepPoint(&aStepPoint);
         aTrack->SetStep(&aStep);
         G4cout << " current energy: " << incomingEnergy
                << " of particle " 
                << aParticle->GetDefinition()->GetParticleName() 
                << " in material " << theMaterials[0]->GetName() << G4endl;
         aParticleChange = 
            (G4ParticleChange*)(theProcesses[i]->PostStepDoIt(*aTrack, aStep));
         G4cout << "NUMBER OF SECONDARIES: "
                << aParticleChange->GetNumberOfSecondaries() << G4endl;
         G4double Tfinal = aParticleChange->GetEnergyChange();
         G4ThreeVector Pfinal = 
                        *(aParticleChange->GetMomentumDirectionChange())
                *sqrt(Tfinal*Tfinal + Tfinal*aParticleChange->GetMassChange());
         G4cout << "FINAL STATE kinetic energy (GeV) = " << Tfinal/GeV 
                << G4endl;
         G4cout << "FINAL STATE momentum (GeV/c) = " 
                << Pfinal*(1./GeV) << G4endl;
         G4Track* second;
         const G4DynamicParticle* aSec;
         G4double eventEnergy;
         G4int Nsec = aParticleChange->GetNumberOfSecondaries();
         HNsec->accumulate((G4float)Nsec);
         if (Nsec == 1) {
            HPx1->accumulate(Pfinal.x()/GeV);
            HPy1->accumulate(Pfinal.y()/GeV);
            HPz1->accumulate(Pfinal.z()/GeV);
            HAng1->accumulate(
                      Pfinal.angle(theDirection)
                      /degree);
            second = aParticleChange->GetSecondary(0);
            HPx2->accumulate(second->GetMomentum().x()/GeV);
            HPy2->accumulate(second->GetMomentum().y()/GeV);
            HPz2->accumulate(second->GetMomentum().z()/GeV);
            HAng2->accumulate(
                      second->GetMomentum().angle(theDirection)
                      /degree);
            eventEnergy = aParticleChange->GetEnergyChange() +
                          aParticleChange->GetMassChange() +
                          second->GetTotalEnergy();
         }
         else if (Nsec == 2) {
            second = aParticleChange->GetSecondary(0);
            HPx1->accumulate(second->GetMomentum().x()/GeV);
            HPy1->accumulate(second->GetMomentum().y()/GeV);
            HPz1->accumulate(second->GetMomentum().z()/GeV);
            HAng1->accumulate(
                   second->GetMomentum().angle(theDirection)
                   /degree);
            eventEnergy = second->GetTotalEnergy();
            second = aParticleChange->GetSecondary(1);
            HPx2->accumulate(second->GetMomentum().x()/GeV);
            HPy2->accumulate(second->GetMomentum().y()/GeV);
            HPz2->accumulate(second->GetMomentum().z()/GeV);
            HAng2->accumulate(
                   second->GetMomentum().angle(theDirection)
                   /degree);
            eventEnergy += second->GetTotalEnergy();
         }

         for (G4int isec = 0; isec < aParticleChange->GetNumberOfSecondaries();
              isec++) {
            second = aParticleChange->GetSecondary(isec);
            aSec = second->GetDynamicParticle();
            G4cout << "SECONDARY " << isec+1 << ": ";
            G4cout << aSec->GetKineticEnergy()/GeV;
            G4cout << aSec->GetMomentum()*(1./GeV);
            G4cout << G4endl;
         }

         G4cout << "TOTAL ENERGY IN EVENT = " << eventEnergy/GeV << G4endl;
         aParticleChange->Clear();

      }  // event loop
   }  // particle loop

   return EXIT_SUCCESS;

}
