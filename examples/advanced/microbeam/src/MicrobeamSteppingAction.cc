// -------------------------------------------------------------------
// $Id: MicrobeamSteppingAction.cc,v 1.3 2006-06-01 22:25:20 sincerti Exp $
// -------------------------------------------------------------------

#include "G4SteppingManager.hh"

#include "MicrobeamSteppingAction.hh"
#include "MicrobeamRunAction.hh"
#include "MicrobeamDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamSteppingAction::MicrobeamSteppingAction(MicrobeamRunAction* run,MicrobeamDetectorConstruction* det)
:Run(run),Detector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamSteppingAction::~MicrobeamSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamSteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 

// COUNT GAS DETECTOR HITS

if (       ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "_CollDet_yoke_")
        &&  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Isobutane")
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha"))
	
        || 
	   ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "_CollDet_gap4_")
        &&  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Isobutane")
	&&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha"))

	||

    	   ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "_CollDet_gap5_")
        &&  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Isobutane")
	&&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha"))

    )
	{
	  Run->AddNbOfHitsGas();	
	}
	
// STOPPING POWER AND BEAM SPOT SIZE AT CELL ENTRANCE

if (       ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Polyprop")
        &&  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "KGM")
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha"))
	
        || 
	   ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Polyprop")
        &&  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalCytoplasm")
	&&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha"))
    )
	{
	
	if(aStep->GetTotalEnergyDeposit()>0) 
	 {
	  FILE *myFile;
	  myFile=fopen("stoppingPower.txt","a");	
	  fprintf(myFile,"%e \n",
	  (aStep->GetTotalEnergyDeposit()/keV)/(aStep->GetStepLength()/micrometer));	
	  fclose (myFile);
	 }

	 G4StepPoint* p1 = aStep->GetPreStepPoint();
         G4ThreeVector coord1 = p1->GetPosition();

         const G4AffineTransform transformation = p1->GetTouchable()->GetHistory()->GetTopTransform();
         G4ThreeVector localPosition = transformation.TransformPoint(coord1);

	 FILE *myFile;
	 myFile=fopen("beamPosition.txt","a");
	 fprintf 
	 ( 
	  myFile,"%e %e \n",
	  localPosition.x()/micrometer,
	  localPosition.y()/micrometer
	 );
	 fclose (myFile);				

	}

// ALPHA RANGE

if (

	(aStep->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "alpha")
	
	&&
	
	(aStep->GetTrack()->GetKineticEnergy()<1e-6)
	
	&&
			
          ( (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalCytoplasm")
	||  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "KGM")	
	||  (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalNucleus") )	
		
   )
	
	{		
	 FILE *myFile;
	 myFile=fopen("range.txt","a");
	 fprintf 
	 ( myFile,"%e %e %e\n",
	  (aStep->GetTrack()->GetPosition().x())/micrometer,
	  (aStep->GetTrack()->GetPosition().y())/micrometer,
	  (aStep->GetTrack()->GetPosition().z())/micrometer
	 );
	 fclose (myFile);				
 	}

// TOTAL DOSE DEPOSIT AND DOSE DEPOSIT WITHIN A PHANTOM VOXEL

if ( 	   (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "physicalNucleus")
  	&& (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalNucleus")
  	&& (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha") )
	{ 
   	 G4double dose = (e_SI*(aStep->GetTotalEnergyDeposit()/eV))/(Run->GetMassNucleus());
   	 Run->AddDoseN(dose);

	 G4ThreeVector v;
    	 Run->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
	}

if ( 	   (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "physicalCytoplasm")
  	&& (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalCytoplasm")
  	&& (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha") )
	{ 
   	 G4double dose = (e_SI*(aStep->GetTotalEnergyDeposit()/eV))/(Run->GetMassCytoplasm());
   	 Run->AddDoseC(dose);

	 G4ThreeVector v;
    	 Run->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
 	}
}
