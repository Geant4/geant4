//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "ExecProcessLevel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"

#include "G4Material.hh"
#include "G4IsotopeVector.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"

#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"

#include "G4HadronInelasticDataSet.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"

#include "G4Gamma.hh"

#include "G4ForceCondition.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"

#include "G4TrackingManager.hh"
#include "G4Region.hh"

#include "TstTarget.hh"
#include "TstDiscreteProcessReader.hh"
// #include "ProcessWrapper.hh"
#include "FTFPWrapper.hh"
#include "QGSPWrapper.hh"


ExecProcessLevel::ExecProcessLevel()
   : ExecBase(),
     fTarget(0), fRegion(0), 
     fProcWrapper(0), fProcManager(0),
     fBeam(0),
     fTrack(0), fStep(0), fStepPoint(0), 
     fPartChange(0)
{
}


ExecProcessLevel::~ExecProcessLevel()
{

   if ( fTarget ) 
   {
      delete fTarget;
      fTarget=0;
   }
   if ( fRegion ) 
   {
      delete fRegion;
      fRegion=0;
   }
   if ( fProcWrapper ) 
   {
      delete fProcWrapper;
      fProcWrapper=0;
   }
   if ( fBeam ) 
   {
      delete fBeam; 
      fBeam=0;
   }
   if ( fTrack ) 
   {
      delete fTrack;
      fTrack=0;
   }
   //if ( fStep ) delete fStep;
   //if ( fStepPoint ) delete fStepPoint;

}

void ExecProcessLevel::InitSetup( const TstReader* pset )
{

   fTarget = new TstTarget();
   
   G4Material* mat = fTarget->ResetMaterial( pset->GetTargetMaterial() );
   
   G4ThreeVector targetSize = pset->GetTargetSize();
   fTarget->SetDimentions( targetSize.x(), targetSize.y(), targetSize.z() );
      
   fTarget->ResetGeom();
   G4VPhysicalVolume* targetPhys = fTarget->ConstructTarget();   
      
   fRegion = new G4Region("Region"); // needed by tracking manager
   targetPhys->GetLogicalVolume()->SetRegion(fRegion);
   fRegion->AddRootLogicalVolume( targetPhys->GetLogicalVolume() );

   return;

}

void ExecProcessLevel::InitBeam( const TstReader* pset )
{

   // fProcWrapper = (dynamic_cast<const TstDiscreteProcessReader*>(pset))->GetProcess();
   
   InitProcess( pset );
   
   assert( fProcWrapper );
   
   assert( fTarget->GetCurrentMaterial() );
   fProcWrapper->InitTarget( fTarget->GetCurrentMaterial() );
      
   if ( pset->GetBeamParticle() == "proton" ) fProcManager = new G4ProcessManager(G4Proton::Proton());
   if ( pset->IsDiscreteProcess() )
   {
      fProcManager->AddDiscreteProcess(fProcWrapper);
   }

   G4ParticleDefinition* partDef = (G4ParticleTable::GetParticleTable())->FindParticle(pset->GetBeamParticle());
   G4double partMass = partDef->GetPDGMass();
   G4double partMom = pset->GetBeamMomentum();
   G4double partEnergy = std::sqrt( partMom*partMom + partMass*partMass );
   
   if ( !fBeam ) fBeam = new Beam();
   fBeam->SetBeam( pset->GetBeamParticle(), partMass, partEnergy ); 
   
   const G4Element* elm = fTarget->GetCurrentMaterial()->GetElement(0);
   G4int A = (G4int)(elm->GetN()+0.5);
   G4int Z = (G4int)(elm->GetZ()+0.5);
   G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);

   fBeam->SetLabV( G4LorentzVector( 0., 0., 
                                    std::sqrt(partEnergy*(partEnergy+2.0*partMass))/GeV,
                                    (partEnergy+partMass+amass)/GeV) );
   
       // if under an agnle, then like this:
       //labv = G4LorentzVector(mom.x()/CLHEP::GeV, mom.y()/CLHEP::GeV, 
       //			      mom.z()/CLHEP::GeV, (e0+mass+amass)/CLHEP::GeV);
       
    fBeam->SetLabP( G4LorentzVector( 0., 0., 
                                     std::sqrt(partEnergy*(partEnergy+2.0*partMass))/GeV,
                                     (partEnergy+partMass+G4Proton::Proton()->GetPDGMass())/GeV) );
         
   // Track
   fTrack = new G4Track( new G4DynamicParticle( partDef, pset->GetDirection(), partEnergy ), 
                         pset->GetTime(), pset->GetPosition() ); 
   fTrack->SetTouchableHandle( G4TouchableHandle( new G4TouchableHistory() ) );

    // Step
    fStep = new G4Step();
    fStep->SetTrack(fTrack);
    fTrack->SetStep(fStep);
    
    // step points
    // pre-step
    fStepPoint = new G4StepPoint();
    fStepPoint->SetPosition( pset->GetPosition() );
    fStepPoint->SetMaterial( fTarget->GetCurrentMaterial() );
    fStepPoint->SetSafety( 10000.*CLHEP::cm );
    fStep->SetPreStepPoint( fStepPoint );
    // post-step
    G4StepPoint* bPoint;
    bPoint = fStepPoint;
    G4ThreeVector bPosition = pset->GetDirection() * pset->GetStep();
    bPosition += pset->GetPosition();
    bPoint->SetPosition(bPosition);
    fStep->SetPostStepPoint( bPoint );
    fStep->SetStepLength( pset->GetStep() );

   return;

}

void ExecProcessLevel::InitProcess( const TstReader* pset )
{

   // G4String name = (dynamic_cast<const TstDiscreteProcessReader*>(pset))->GetProcessName();
   G4String name = (pset)->GetPhysics();

   if ( name.find("ftf") != std::string::npos )
   {
      fProcWrapper = new FTFPWrapper();
   }
   else if ( name.find("qgs") != std::string::npos )
   {
       fProcWrapper = new QGSPWrapper();
   }
    
   if (!fProcWrapper) 
   { 
      G4cout 
	     << " generator " << name << " is unavailable"
	     << G4endl;
      exit(1);
   } 
    
   std::cout << " process = " << fProcWrapper->GetProcessName() << std::endl;
    
   if ( name.find("lund-str-fragm") != std::string::npos )
   {
      fProcWrapper->UseG4LundStringFragm(true);
   }
    
   fProcWrapper->Compose();
   
   return;

}

G4VParticleChange* ExecProcessLevel::DoEvent()
{

   fTrack->SetKineticEnergy( fBeam->GetBeamEnergy()-fBeam->GetBeamPartMass() );
   
   fPartChange = fProcWrapper->PostStepDoIt( *fTrack, *fStep );
   
   return fPartChange;

}
