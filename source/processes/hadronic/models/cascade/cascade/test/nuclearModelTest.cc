#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"
#include <string.h>

int main(){

   G4NucleusModel carbonNucleus;
   carbonNucleus.CreateModel(12,6);
  
  cout << "Carbon nucleus created" << endl;

  cout << "Number of protons: " << carbonNucleus.GetZ() << endl;
  cout << "Number of nucleons: " << carbonNucleus.GetA() << endl;
  cout << "Mass of carbon nucleus: " << carbonNucleus.GetNuclearMass()<< " MeV" << endl;
  cout << "Atomic mass by Weitzsaecker's mass formula: " << carbonNucleus.GetAtomicMass() <<" MeV" <<  endl;
  cout << "Radius of nucleus: " << carbonNucleus.GetRadius() << " fermi" << endl;
  cout << "Proton density of the nucleus: " << carbonNucleus.GetProtonDensity(0) << " 1/m3 " << endl;
  cout << "Neutron density of the nucleus: " << carbonNucleus.GetNeutronDensity(0) <<" 1/m3" << endl;
  cout << "Initial excitation energy: " << carbonNucleus.GetExcitationEnergy() << " MeV" << endl;
 
  cout << "Let's excite the nucleus!" << endl;
  cout << "-----------------------------" << endl;
  carbonNucleus.ChangeExcitationEnergy(10);
  cout << "Excitation energy now: " << carbonNucleus.GetExcitationEnergy() <<" MeV" <<  endl;

  cout <<"Let's collide the nucleus!" << endl;
  cout << "-----------------------------" << endl;
  G4DynamicParticle* incidentParticle = new G4DynamicParticle;
  incidentParticle->SetDefinition(G4Proton::Proton());
  G4ThreeVector incidentMomentum(100, 100, 100);
  G4double protonMass = G4Proton::Proton()->GetPDGMass();
  G4double totalEnergy = sqrt(sqr(incidentMomentum.mag()) + sqr(protonMass));
  G4LorentzVector incident4Momentum(100.0, 100.0, 100.0, totalEnergy);
  incidentParticle->Set4Momentum(incident4Momentum);

  carbonNucleus.DoCollision(incidentParticle);
  G4ThreeVector interactionPoint = carbonNucleus.GetInteractionPoint();
 

  if(carbonNucleus.IsInside(interactionPoint)){
    cout << "Interaction happened!" << endl;

    int i;

    cout << "interaction point: " << interactionPoint <<" fermi" <<  endl;

    G4DynamicParticle *targetParticle = carbonNucleus.ReturnTargetNucleon();
    cout << "Struck particle was: "<< targetParticle->GetDefinition()->GetParticleName() << endl;
    
    G4ThreeVector p = targetParticle->GetMomentum();
      cout << "Struck particle's momentum: " << p  <<  endl;
  }
  else
    cout << "No interaction!" << endl;
  
  carbonNucleus.PrintModel();
  return 0;
}



