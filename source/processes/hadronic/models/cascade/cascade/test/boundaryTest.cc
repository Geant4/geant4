
#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"

//Test for reflections and transitions between regions.
//The private functions in G4NucleusModel must be public during the test.

int main(){
  G4NucleusModel* nucleus = new G4NucleusModel();
  G4VRegionModel* firstRegion;
  G4DynamicParticle* incidentParticle = new G4DynamicParticle();
  G4DynamicParticle* incidentParticle2 = new G4DynamicParticle();
  G4DynamicParticle* incidentParticle3 = new G4DynamicParticle();
  G4DynamicParticle* incidentParticle4 = new G4DynamicParticle();
  G4DynamicParticle* incidentParticle5 = new G4DynamicParticle();
  G4ThreeVector initialPoint(0.0, 0.0, 0.0);
  G4bool crossed;
  G4int i;

  //testing CreateModel / PrintModel

  nucleus->CreateModel(18, 8);
  nucleus->PrintModel();
  G4double radius1 = nucleus->GetRegionRadius(1);
  G4double radius2 = nucleus->GetRegionRadius(2);
  G4double radius3 = nucleus->GetRegionRadius(3);  

  firstRegion = nucleus->GetRegion(0);
  incidentParticle->SetDefinition(G4Proton::Proton());
  incidentParticle2->SetDefinition(G4Proton::Proton());
  incidentParticle3->SetDefinition(G4Proton::Proton());
  incidentParticle4->SetDefinition(G4Proton::Proton());
  incidentParticle5->SetDefinition(G4Proton::Proton());
  G4double protonMass = G4Proton::Proton()->GetPDGMass();
  G4ThreeVector momentum(1.0, 0.0, 0.0);
  incidentParticle->SetMomentum(momentum);
  nucleus->DoCollision(incidentParticle);

  cout << "Nucleus collided" << endl;
  G4double radius = nucleus->GetRadius();

  //testing BoundaryCrossed

  crossed = nucleus->BoundaryCrossed(firstRegion, initialPoint, incidentParticle, radius3/2);
  if(crossed) cout << "Boundary crossed (region1 -> region2)!" << endl;
  else cout << "No boundary crossed!" << endl;
  cout << endl;

  G4ThreeVector point2(0.0, 0.0, 0.5*radius2+0.5*radius3);
  G4VRegionModel* thirdRegion = nucleus->GetRegion(point2.mag());
  G4ThreeVector momentum2(0.0, 0.0, -1.0);
 
  incidentParticle2->SetMomentum(momentum2);
  nucleus->DoCollision(incidentParticle2);
  crossed = nucleus->BoundaryCrossed(thirdRegion, point2, incidentParticle2, 2*radius3/3);
  if(crossed) cout << "Boundary crossed (region3 -> region2)!" << endl << endl;
  else cout << "No boundary crossed!" << endl;
  
  G4ThreeVector point3(0.0, 0.0, 0.5*radius1 + 0.5*radius2);
  G4VRegionModel* secondRegion = nucleus->GetRegion(point3.mag());
  G4ThreeVector momentum3(0.0, 0.0, -1.0);
   incidentParticle3->SetMomentum(momentum3);
  nucleus->DoCollision(incidentParticle3);

  crossed = nucleus->BoundaryCrossed(secondRegion, point3, incidentParticle3, 0.4*radius3);
  if(crossed) cout << "Boundary crossed (region2 -> region1)!" << endl << endl;
  else cout << "No boundary crossed!" << endl;

  
  G4ThreeVector momentum4(0.0, 0.0, 1.0);
  incidentParticle4->SetMomentum(momentum4);
  nucleus->DoCollision(incidentParticle);

  crossed = nucleus->BoundaryCrossed(secondRegion, point3, incidentParticle4, 0.4*radius3);
  if(crossed) cout << "Boundary crossed (region2 -> region3)!" << endl << endl;
  else cout << "No boundary crossed!" << endl;

  crossed = nucleus->BoundaryCrossed(thirdRegion, point2, incidentParticle4, 0.4*radius3);
  if(crossed) cout << "Boundary crossed (region3 -> region4)!" << endl << endl;
  else cout << "No boundary crossed!" << endl;
  
  cout << "------------------------------------------------" << endl;

  //testing GetNextRegion

  G4VRegionModel* nextRegion;

  nextRegion = nucleus->GetNextRegion(firstRegion, initialPoint, incidentParticle);
  cout << "next Region 1 radius: " << nextRegion->GetOuterRadius() << endl << endl;

  nextRegion = nucleus->GetNextRegion(thirdRegion, point2, incidentParticle2);
  cout << "next Region 2 radius: " << nextRegion->GetOuterRadius() << endl << endl;
 
  nextRegion = nucleus->GetNextRegion(secondRegion, point3, incidentParticle3);
  cout << "next Region 3 radius: " << nextRegion->GetOuterRadius() << endl << endl;
 
  nextRegion = nucleus->GetNextRegion(secondRegion, point3, incidentParticle4);
  cout << "next Region 4 radius: " << nextRegion->GetOuterRadius() << endl << endl; 
  
  cout << "--------------------------------------------------" << endl;
  //testing GetCrossingPoint
  G4ThreeVector crossingPoint1 = nucleus->GetCrossingPoint(firstRegion, initialPoint, incidentParticle);
  cout << "crossing point 1: " << crossingPoint1 << endl << endl;

  G4ThreeVector crossingPoint2 = nucleus->GetCrossingPoint(thirdRegion, point2, incidentParticle2);
  cout << "crossing point 2: " << crossingPoint2 << endl << endl;

  G4ThreeVector crossingPoint3 = nucleus->GetCrossingPoint(secondRegion, point3, incidentParticle3);
  cout << "crossing point 3: " << crossingPoint3 << endl << endl;

  G4ThreeVector crossingPoint4 = nucleus->GetCrossingPoint(secondRegion, point3, incidentParticle4);
  cout << "crossing point 4: " << crossingPoint4 << endl << endl;

  G4ThreeVector momentum5(1.0, 1.0, 1.0);
  incidentParticle5->SetMomentum(momentum5);
  G4ThreeVector crossingPoint5 = nucleus->GetCrossingPoint(firstRegion, initialPoint, incidentParticle5);
  cout << "crossing point 5: " << crossingPoint5 << endl;
  cout << "should have been: " << firstRegion->GetOuterRadius()/sqrt(3) << endl << endl; 

  cout << "---------------------------------------------------" << endl;

  //testing DoBoundaryTransition
  
  momentum[0] = 1000.0;
  G4double energy1 = sqrt(sqr(momentum.mag()) + sqr(protonMass));
  cout << "momentum.mag(): " << momentum.mag() << endl;
  cout << "energy1: " << energy1 << endl;
  G4LorentzVector fourMomentum1(1000.0, 0.0, 0.0, energy1);
  incidentParticle->Set4Momentum(fourMomentum1);
  cout << "initial momentum 1: " << incidentParticle->Get4Momentum() << endl;
  nextRegion = nucleus->DoBoundaryTransition(firstRegion, crossingPoint1, incidentParticle);
  cout << "Region radius after boundary transition 1: " << nextRegion->GetOuterRadius() << endl; 
  cout << "final momentum 1: " << incidentParticle->Get4Momentum() << endl << endl;

  momentum2 = (0.0, 0.0, -100.0);
  G4double energy2 = sqrt(sqr(momentum2.mag()) + sqr(protonMass));
  cout << "momentum2.mag(): " << momentum2.mag() << endl;
  cout << "energy2: " << energy2 << endl;
  G4LorentzVector fourMomentum2(0.0, 0.0, -100.0, energy2);
  incidentParticle2->Set4Momentum(fourMomentum2);
  cout << "initial momentum 2: " << incidentParticle2->Get4Momentum() << endl;
  nextRegion = nucleus->DoBoundaryTransition(thirdRegion, crossingPoint2, incidentParticle2);
  cout << "Region radius after boundary transition 2: " << nextRegion->GetOuterRadius() << endl; 
  cout << "final momentum 2: " << incidentParticle2->Get4Momentum() << endl << endl;  
  

  momentum5 = (30.0, 30.0, 30.0);
  G4double energy5 = sqrt(sqr(momentum5.mag()) + sqr(protonMass));
  cout << "momentum5.mag(): " << momentum5.mag() << endl;
  cout << "energy5: " << energy5 << endl;
  G4LorentzVector fourMomentum5(30.0, 30.0, 30.0, energy5);
  incidentParticle5->Set4Momentum(fourMomentum5);
  cout << "initial momentum 5: " << incidentParticle5->Get4Momentum() << endl;
  nextRegion = nucleus->DoBoundaryTransition(firstRegion, crossingPoint5, incidentParticle5);
  cout << "Region radius after boundary transition 5: " << nextRegion->GetOuterRadius() << endl; 
  cout << "final momentum 5: " << incidentParticle5->Get4Momentum() << endl << endl;  
  cout << "-----------------------------------------------------------" << endl;

  //testing CalculateInteractionPoint
  
  
  G4ThreeVector mom(1.0, 0.0, -1.0);
  G4double energy = sqrt(sqr(mom.mag()) + sqr(protonMass));
  G4LorentzVector fourMom(1.0, 0.0, -1.0, energy);
  incidentParticle->Set4Momentum(fourMom);
  nucleus->DoCollision(incidentParticle);
  G4double lambda = nucleus->GetInteractionLength();
  cout << "lambda: " << lambda << endl;
  G4ThreeVector interactionPoint1 = nucleus->CalculateInteractionPoint();
  cout << "interaction point 1: " << interactionPoint1 << endl << endl;

  nucleus->DoCollision(incidentParticle2);
  G4ThreeVector interactionPoint2 = nucleus->CalculateInteractionPoint();
  cout << "interaction point 2: " << interactionPoint2 << endl;

  //testing DoCollision
  //nucleus->DoCollision(incidentParticle);
  return 0;
}









