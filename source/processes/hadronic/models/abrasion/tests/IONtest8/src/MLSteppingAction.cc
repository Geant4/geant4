////////////////////////////////////////////////////////////////////////////////
//
#include "MLSteppingAction.hh"

#include "G4ios.hh"
#include <assert.h>

#include "MLSD.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"
#include "G4StepStatus.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "MLGeometryConstruction.hh"
#include "MLAnalysisManager.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLSteppingAction::MLSteppingAction (MLGeometryConstruction* det)
  :geometry(det)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLSteppingAction::~MLSteppingAction ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MLSteppingAction::UserSteppingAction (const G4Step* fStep)
{
  if (fStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    const G4VTouchable* preStepTouchable =
      fStep->GetPreStepPoint()->GetTouchable();
    const G4VTouchable* postStepTouchable =
      fStep->GetPostStepPoint()->GetTouchable();

    if (postStepTouchable->GetVolume()) {
      G4int layer1 = 0;
      if (preStepTouchable->GetVolume()->GetName() != "World") {
        layer1 = preStepTouchable->GetVolume()->GetCopyNo();
      } else {
        layer1 = -1;
      }
      G4int layer2 = 0;
      if (postStepTouchable->GetVolume()->GetName() != "World") {
        layer2 = postStepTouchable->GetVolume()->GetCopyNo();
      } else {
        layer2 = -1;
      }

      G4bool isFDet1    = geometry->IsAFluxDetector(layer1);
      G4bool isFDet2    = geometry->IsAFluxDetector(layer2);

      G4int layerNumber = -9999;
      if (layer1 != -1 && layer2 != -1) {
        if (layer2 > layer1 && isFDet2) {
          layerNumber = layer2;                  //In geometry propagating forward into flux detector boundary
        } else if (layer2 < layer1 && isFDet1) {
          layerNumber = layer1;                  //In geometry propagating backwards into flux detector boundary
        }
      } else if (layer1 == -1 && isFDet2) {
        layerNumber = 0;                         //Outside geometry propaging into flux detector baoundary - this is the source
      } else if (layer2 == -1) {
        if (layer1 == 0 && isFDet1) {
          layerNumber = 0;                       //Inside geometry propagating back through front surface
        } else if (layer1 == geometry->GetNbOfLayers()-1 &&
		   geometry->IsAFluxDetector(layer1+1)) {
          layerNumber = layer1+1;                  //Inside geometry propagating out of the back of the shield.
	}
      }

      if (layerNumber >=0) {

        G4Track* fTrack = fStep->GetTrack();
        G4int pName = fTrack->GetDefinition()->GetPDGEncoding();
        if ( pName == 2212 || pName == 2112 || // proton and neutron
             pName == 11   || pName == 22   || // electron & photon
             pName == 13   || pName == -13  || // muon- & muon+
             pName == 211  || pName == -211) { // poin+ & pion-
          G4double energy           = fTrack->GetKineticEnergy();
          G4ThreeVector direction   = fTrack->GetMomentumDirection() ;
          G4ThreeVector position    = fTrack->GetPosition() ;
          G4ThreeVector translation = postStepTouchable->GetVolume()
                                         ->GetObjectTranslation();
          position -= translation;
	    //	    const G4VTouchable* preStepTouchable = fStep->GetPostStepPoint()->GetTouchable();
          G4ThreeVector normal(0.,0.,0.);
          if (layer2 >= 0) {
            normal = postStepTouchable->GetVolume()->GetLogicalVolume()
              ->GetSolid()->SurfaceNormal(position);
          } else {
            normal = -preStepTouchable->GetVolume()->GetLogicalVolume()
              ->GetSolid()->SurfaceNormal(position);
          }
          G4double theta = direction.angle(normal);
	    //	    if (theta < 90. && pName == 2112) G4cout << postStepTouchable->GetVolume()->GetName() << normal << position << theta << endl;
          if ((theta > 90.*deg && fStep->GetStepLength() != 0.) &&// Make sure the particle is indeed entering the volume!
              (geometry->GetShape()==SPHERE ||
                (geometry->GetShape()==SLAB &&
                  (position.x() < geometry->GetWorldSizeXY() &&
                   position.x() >-geometry->GetWorldSizeXY() &&
                   position.y() < geometry->GetWorldSizeXY() &&
                   position.y() >-geometry->GetWorldSizeXY())))) {
            if (geometry->GetShape() == SLAB) {
              theta = direction.getTheta();
            } else if (position.mag() < geometry->GetLayerRadius(layer2)) {
              theta = 180.*deg - theta;
            }

            G4ThreeVector aparticle = G4ThreeVector(pName,energy,theta);
            G4double weight         = fTrack->GetWeight();
            G4SDManager* SDman      = G4SDManager::GetSDMpointer();
            MLSD* sd                = dynamic_cast<MLSD*>
              (SDman->FindSensitiveDetector("detectorSD") );
            MLHit* aHit             = new MLHit();
            aHit->SetParticle(layerNumber,weight,aparticle);
            sd->MLCollection->insert(aHit);
          }
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
