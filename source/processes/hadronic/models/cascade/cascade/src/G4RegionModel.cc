#include "../include/G4RegionModel.hh"


//initializes static variable numberOfRegions
G4int G4RegionModel::numberOfRegions = 0;

G4RegionModel::G4RegionModel(){ 
   numberOfRegions++; 
   myRegionNumber = numberOfRegions;
}  

G4RegionModel::~G4RegionModel(){}

void G4RegionModel::InitializeRadius(){ 
  //radiuses calculated from density formula used in InitializeDensity

  const G4double regionDensities[3] = {0.9, 0.2, 0.1}; //proportions of the center density value 
  const G4double a = 0.6;
  G4double b, c; //help variables

 if(numberOfRegions == 1){
   b = 1 / regionDensities[2];
   innerRadius = 0.0;
  }
    
 if(numberOfRegions == 3){
    b = 1/regionDensities[myRegionNumber-1];
    if(myRegionNumber >= 2) c = 1/regionDensities[myRegionNumber-2];
 }
  outerRadius = a*log((b - 1.0)*exp(nuclearRadius/a) + b);
  if(myRegionNumber >= 2) innerRadius = a*log((c - 1.0)*exp(nuclearRadius/a) + c);
  else innerRadius = 0.0;
}		   


void G4RegionModel::InitializeDensity(){
  // densities based on density formula from P.E.Hodgson;
  // Introductory Nuclear Physics

  const G4double a = 0.6;
  const G4double rho0 = 0.08;
  const G4double protonRatio = (nucleusZ + 0.0) / nucleusA;
  const G4double neutronRatio = (nucleusA - nucleusZ + 0.0) / nucleusA;
  G4double matterDensity;

    matterDensity = rho0 / (1 + exp(((outerRadius + innerRadius)/2 - nuclearRadius)/a));
    protonDensity = protonRatio*matterDensity;
    neutronDensity = neutronRatio*matterDensity;
}
    
void G4RegionModel::InitializeFermi(){
  //initializes fermi momenta and energies
  G4double protonMass = G4Proton::Proton()->GetPDGMass();
  G4double neutronMass = G4Neutron::Neutron()->GetPDGMass();

  protonFermiMomentum = GetFermiMomentum(protonDensity, protonMass);
  neutronFermiMomentum = GetFermiMomentum(neutronDensity, neutronMass);
  
  protonFermiEnergy = sqr(protonFermiMomentum)/(2*protonMass);
  neutronFermiEnergy = sqr(neutronFermiMomentum)/(2*neutronMass);
}

void G4RegionModel::InitializePotentialEnergy(){
  const G4double BE = 7.0*MeV; //binding energy of the most loosely bound nucleon
  protonPotentialEnergy = -protonFermiEnergy - BE;
  neutronPotentialEnergy = -neutronFermiEnergy - BE;
}

void G4RegionModel::CreateRegion(G4int A, G4int Z){
  
  const G4double radius0 = 1.2;
  const G4double oneThird = 1./3.;
 
  nucleusA = A; 
  nucleusZ = Z;
  nuclearRadius =  radius0*pow(nucleusA, oneThird);

  InitializeRadius();
  InitializeDensity();
  InitializeFermi();
  InitializePotentialEnergy();
} 

G4double G4RegionModel::GetFermiMomentum(G4double density,  
					 G4double mass){
    const G4double oneThird = 1./3.;
    const G4double Hbarc =197.32696; //MeV*fm, from PDG 

  return (Hbarc*pow(3*pi2*density, oneThird));
  //formula from class G4FermiMomentum
}


G4double G4RegionModel::GetOuterRadius(){
  return outerRadius;
}

G4double G4RegionModel::GetProtonDensity(){
  return protonDensity;
}

G4double G4RegionModel::GetNeutronDensity(){
  return neutronDensity;
}

G4double G4RegionModel::GetProtonPotentialEnergy(){
  return protonPotentialEnergy;
}

G4double G4RegionModel::GetNeutronPotentialEnergy(){
  return neutronPotentialEnergy;
}

G4double G4RegionModel::GetProtonMaximumMomentum(){
  return protonFermiMomentum; 
}

G4double G4RegionModel::GetProtonMaximumEnergy(){
  return protonFermiEnergy;
}

G4double G4RegionModel::GetNeutronMaximumEnergy(){
  return neutronFermiEnergy;
}

G4double G4RegionModel::GetNeutronMaximumMomentum(){
  return neutronFermiMomentum;
}

void G4RegionModel::PrintRegion(){ 
  cout << "Region number: " << myRegionNumber  << endl; 
  cout << "Outer radius: " << outerRadius <<" fm" <<  endl; 
  cout << "Proton density: " << protonDensity << " protons/fm3" << endl; 
  cout << "Neutron density: " << neutronDensity << " neutrons/fm3" << endl;  
  cout << "Proton Fermi momentum: " << protonFermiMomentum <<" MeV" << endl;
  cout << "Neutron Fermi momentum: " << neutronFermiMomentum << " MeV" << endl;  
  cout << "Proton Fermi energy: " << protonFermiEnergy <<" MeV" << endl;
  cout << "Neutron Fermi energy: " << neutronFermiEnergy <<" MeV" << endl;    
  cout << "Proton potential energy: " << protonPotentialEnergy <<" MeV" << endl; 
  cout << "Neutron potential energy: " << neutronPotentialEnergy <<" MeV" << endl; 
}

G4DynamicParticle* G4RegionModel::GenerateProton(){
    
    G4DynamicParticle *particle = new G4DynamicParticle();
    particle->SetDefinition(G4Proton::Proton());
  
    G4double cosMu  = 1 - 2 * G4UniformRand(); 
    G4double phi = 2 * pi * G4UniformRand(); //azimuthal angle
    G4double mu = acos(cosMu); //polar angle

    G4double q = protonFermiMomentum * sqrt(G4UniformRand());  //magnitude of the momentum
    G4ThreeVector p(q*cos(phi)*sin(mu), q*sin(phi)*sin(mu), q*cos(mu));  //spherical coordinates 
    G4double totalEnergy = sqrt(sqr(p.mag()) + sqr(particle->GetMass()));
    G4LorentzVector fourMomentum(q*cos(phi)*sin(mu), q*sin(phi)*sin(mu), q*cos(mu), totalEnergy);
    
    particle->SetMomentum(p);
    particle->Set4Momentum(fourMomentum); 
    
  return particle;
}

G4DynamicParticle* G4RegionModel::GenerateNeutron(){
    G4DynamicParticle *particle = new G4DynamicParticle();
    particle->SetDefinition(G4Neutron::Neutron());
    
    G4double cosMu  = 1 - 2 * G4UniformRand(); 
    G4double phi = 2 * pi * G4UniformRand(); 
    G4double mu = acos(cosMu); 

    G4double q = neutronFermiMomentum * sqrt(G4UniformRand()); 
    G4ThreeVector p(q*cos(phi)*sin(mu), q*sin(phi)*sin(mu), q*cos(mu)); 
    G4double totalEnergy = sqrt(sqr(p.mag()) + sqr(particle->GetMass()));
    G4LorentzVector fourMomentum(q*cos(phi)*sin(mu), q*sin(phi)*sin(mu), q*cos(mu), totalEnergy);
    
    particle->SetMomentum(p);
    particle->Set4Momentum(fourMomentum); 
    
  return particle;
}











