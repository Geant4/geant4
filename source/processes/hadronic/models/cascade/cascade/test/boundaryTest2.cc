#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"

//Test for reflections and transitions between regions.
//The private functions in G4NucleusModel must be made public for the test.

int main(){
  G4NucleusModel* nucleus = new G4NucleusModel();
  G4DynamicParticle* incidentParticle = new G4DynamicParticle();
  G4double protonMass =  G4Proton::Proton()->GetPDGMass();
  G4VRegionModel* nextZone;

  nucleus->CreateModel(22, 10);
  nucleus->PrintModel();

  G4double radius1 = nucleus->GetRegionRadius(1);
  //  cout << "radius1: " << radius1 << endl;
  G4double radius2 = nucleus->GetRegionRadius(2);
  // cout << "radius2: " << radius2 << endl;
  G4ThreeVector initialPoint(1./2.*radius1 + 1./2.*radius2,
			     1./3.*radius1, 0.0);
  //cout << "initial point: " << initialPoint << endl;
  G4VRegionModel* secondRegion = nucleus->GetRegion(initialPoint.mag());
  //cout << "second region radius: " << secondRegion->GetOuterRadius() << endl;
  G4ThreeVector initialMomentum(-200.0, 0.0, 0.0);
  G4double totalEnergy = sqrt(sqr(initialMomentum.mag()) + sqr(protonMass));
  G4LorentzVector initial4Momentum(-200.0, 0.0, 0.0, totalEnergy);

  incidentParticle->SetDefinition(G4Proton::Proton());
  incidentParticle->Set4Momentum(initial4Momentum);

  //let's do the transitions!
  G4ThreeVector crossingPoint = nucleus->GetCrossingPoint(secondRegion, initialPoint, incidentParticle);
  cout << "initial momentum: " << initial4Momentum << endl << endl;
  nextZone = secondRegion;
  G4VRegionModel* test;

  for(int i=1; i<4; i++){
    cout << "crossing point " << i <<": " << crossingPoint << endl;
    //    test = nucleus->GetNextRegion(nextZone, crossingPoint, incidentParticle);
    // if(test!=NULL) cout <<"next region radius: " <<test->GetOuterRadius() << endl;
    if(nextZone!=NULL) nextZone = nucleus->DoBoundaryTransition(nextZone, crossingPoint, incidentParticle);
    cout << "momentum after " << i << ". boundary: " << incidentParticle->Get4Momentum() << endl << endl; 
    //    cout << "radius after: " << nextZone->GetOuterRadius() << endl;
    if(nextZone!=NULL) crossingPoint = nucleus->GetCrossingPoint(nextZone, crossingPoint, incidentParticle);
  }

  if(nextZone == NULL)
    cout << "Particle out of nucleus!" << endl;

  return 0;
}









