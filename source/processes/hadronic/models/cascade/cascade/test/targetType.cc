#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"
#include <string.h>

//This test tests the type of the generated target particle of class
//G4NucleusModel. The type is tested with 3 different energies numberOfTimes
//times.

int main(){
  G4double E1, E2, E3, *energyTable, totalEnergy;
  ifstream in;
  G4int numberOfTargets, numberOfEnergies, numberOfProtons = 0, numberOfNeutrons = 0, i, j;
  G4ThreeVector interactionPoint;
  G4LorentzVector fourMomentum;
  G4NucleusModel nucleus;
  G4DynamicParticle *targetParticle; 
  G4String particleName;
  G4double protonMass =  G4Proton::Proton()->GetPDGMass();

  nucleus.CreateModel(20,15);
  G4DynamicParticle* incidentParticle = new G4DynamicParticle;
  incidentParticle->SetDefinition(G4Proton::Proton());
  in.open("targetType.in", ios::in); 

  in >> numberOfEnergies;
  energyTable = new G4double[numberOfEnergies];

  for(j=0; j<numberOfEnergies; j++)  in >> energyTable[j];
  
  //cout << "Here we are 2" << endl;
  in >> numberOfTargets;

  for(j=0; j<numberOfEnergies; j++){
    incidentParticle->SetKineticEnergy(energyTable[j]);
    totalEnergy = sqrt(2*energyTable[j]*protonMass + sqr(protonMass));
    fourMomentum[0] = sqrt(2*energyTable[j]*protonMass); 
    fourMomentum[1] = 0; fourMomentum[2] = 0; fourMomentum[3] = totalEnergy;
    incidentParticle->Set4Momentum(fourMomentum);
    //   cout << "fourMomentum: " << fourMomentum << endl; 
    //cout << "hete we are 3" << endl;
    //cout << "number of targets " << numberOfTargets << endl;
    for(int i=0; i<numberOfTargets; i++){
      nucleus.DoCollision(incidentParticle);
      //cout << "i: " << i << " here we are 4 " << endl;
      interactionPoint = nucleus.GetInteractionPoint();
      //cout << "here we are 5 " << endl;
      //cout <<"i: " << i  <<" interaction point: " << interactionPoint << endl;
      if(nucleus.IsInside(interactionPoint)){
	targetParticle = nucleus.ReturnTargetNucleon();
	particleName = targetParticle->GetDefinition()->GetParticleName();
	if(strcmp(particleName, "proton")==0) numberOfProtons++;
	else numberOfNeutrons++;
	//cout << "i: " << i <<  " Here we are 6!" << endl;
      }
    }
    
    cout << "energy: " << energyTable[j] << ", protons: " <<numberOfProtons << ", neutrons: " << numberOfNeutrons << endl;

    numberOfProtons = 0;
    numberOfNeutrons = 0;

  }
  return 0;    
}  
