#include "../include/G4NucleusModel.hh"

//constants
const G4double radius0 = 1.2;
const G4double oneThird = 1./3.;

//initialization of the static pointer
G4BertiniData* G4NucleusModel::xSecDataBase = NULL;

G4NucleusModel::G4NucleusModel(){}

G4NucleusModel::~G4NucleusModel(){
  if(regions.size() > 0){
    for(regionIterator iter = regions.begin(); iter < regions.end(); iter++)  delete(*iter);
  }
  if(xSecDataBase!=NULL) delete xSecDataBase;
}


void G4NucleusModel::SetParameters(const G4int A, const G4int Z){
  myA = A;
  myZ = Z;
  bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z,A);
  radius = radius0*pow(myA, oneThird); 
  atomicMass = CalculateAtomicMass();
}

void G4NucleusModel::SetMomentum(G4ThreeVector p){
  momentum = p;
  kineticEnergy = sqrt(sqr(momentum.mag())+sqr(atomicMass)) - atomicMass;
}

void G4NucleusModel::ChangeParameters(const G4int deltaA, const G4int deltaZ){
  myA += deltaA;
  myZ += deltaZ;
  bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(myZ,myA);
  radius = radius0*pow(myA, oneThird); 
  atomicMass = CalculateAtomicMass();  
}

void G4NucleusModel::CreateModel(const G4int A, const G4int Z){ 
   myA = A; 
   myZ = Z; 
   excitationEnergy = 0.0;
   momentum = (0.0, 0.0, 0.0);
   kineticEnergy = 0.0;
   bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z,A); 
   radius = radius0*pow(myA, oneThird); 
   atomicMass = CalculateAtomicMass();

   for(int i=0; i<NUMBER_OF_REGIONS; i++) 
     regions.push_back(new G4RegionModel());

   regionIterator iter;
   for(iter = regions.begin(); iter < regions.end(); iter++) (*iter)->CreateRegion(A,Z);
}


void G4NucleusModel::ChangeExcitationEnergy(G4double deltaE){
  excitationEnergy += deltaE;
}

G4double G4NucleusModel::GetA(){
  return myA;
}

G4double G4NucleusModel::GetZ(){
  return myZ;
}

G4double G4NucleusModel::GetRadius(){
  return radius;
}

G4double G4NucleusModel::GetRegionRadius(G4int regionNumber){
  if(regionNumber > regions.size() || regionNumber < 1) return 0;
  else{
    regionIterator iter;
    G4int counter = 1;
    for(iter = regions.begin(); iter < regions.end(); iter++){
      if(counter == regionNumber) return (*iter)->GetOuterRadius();
      else counter++;
     }
  }
}

G4double G4NucleusModel::GetExcitationEnergy(){
  return excitationEnergy;
}

G4double G4NucleusModel::GetNuclearMass(){
   return myZ*G4Proton::Proton()->GetPDGMass() + 
          (myA-myZ)*G4Neutron::Neutron()->GetPDGMass() -
          bindingEnergy;
   //mass = sum of proton masses + sum of neutron masses - binding energy
}

G4double G4NucleusModel::GetAtomicMass(){
  return atomicMass;
}

G4double G4NucleusModel::GetBindingEnergy(){
  return bindingEnergy;
}

G4double G4NucleusModel::CalculateAtomicMass(){

    const G4double electron_mass = G4Electron::Electron()->GetPDGMass()/MeV;
    const G4double proton_mass = G4Proton::Proton()->GetPDGMass()/MeV;
    const G4double neutron_mass = G4Neutron::Neutron()->GetPDGMass()/MeV;
    const G4double deuteron_mass = G4Deuteron::Deuteron()->GetPDGMass()/MeV;
    const G4double alpha_mass = G4Alpha::Alpha()->GetPDGMass()/MeV;

    G4int theZ = G4int(myZ + 0.5);
    G4int theA = G4int(myA + 0.5);
      
    if( theA == 1 )
    {
      if( theZ == 0 ) return neutron_mass*MeV;
      if( theZ == 1 ) return proton_mass*MeV + electron_mass*MeV;   // hydrogen
    }
    else if( theA == 2 && theZ == 1 )
    {
      return deuteron_mass*MeV;
    }
    else if( theA == 4 && theZ == 2 )
    {
      return alpha_mass*MeV;
    }
 
    // Weitzsaecker's Mass formula
       
    G4double mass =
      (myA-myZ)*neutron_mass + myZ*proton_mass + myZ*electron_mass
      - 15.67*myA                                          // nuclear volume
      + 17.23*pow(myA,2./3.)                               // surface energy
      + 93.15*pow(myA/2.-myZ,2.)/myA                       // asymmetry
      + 0.6984523*pow(myZ,2.)*pow(myA,-oneThird);          // coulomb
    G4int ipp = (theA - theZ)%2;                           // pairing
    G4int izz = theZ%2;
    if( ipp == izz )mass += (ipp+izz-1) * 12.0 * pow(myA,-0.5);
    return mass*MeV;
}


G4double G4NucleusModel::GetProtonDensity(G4ThreeVector position){
  G4VRegionModel* thisRegion = GetRegion(position.mag());
  return thisRegion->GetProtonDensity();
}

G4double G4NucleusModel::GetNeutronDensity(G4ThreeVector position){
  G4VRegionModel* thisRegion = GetRegion(position.mag());
  return thisRegion->GetNeutronDensity();
}

G4double G4NucleusModel::GetKineticEnergy(){
  return kineticEnergy;
}

G4ThreeVector G4NucleusModel::GetMomentum(){
  return momentum;
}


G4double G4NucleusModel::GetProtonFermiMomentum(G4ThreeVector position){
  G4VRegionModel* thisRegion=GetRegion(position.mag());
  return thisRegion->GetProtonMaximumMomentum();
}

G4double G4NucleusModel::GetNeutronFermiMomentum(G4ThreeVector position){
  G4VRegionModel* thisRegion=GetRegion(position.mag());
  return thisRegion->GetNeutronMaximumMomentum();
}

G4double G4NucleusModel::GetProtonPotentialEnergy(G4ThreeVector position){
  G4VRegionModel* thisRegion=GetRegion(position.mag());
  return thisRegion->GetProtonPotentialEnergy();
}

G4double G4NucleusModel::GetNeutronPotentialEnergy(G4ThreeVector position){
  G4VRegionModel* thisRegion=GetRegion(position.mag());
  return thisRegion->GetNeutronPotentialEnergy();
}


G4ThreeVector G4NucleusModel::GetIncidencePoint(){
  //algorithm from Methods in Subnuclear Physics, Notes on Monte Carlo method 
  G4double b              = radius * sqrt(G4UniformRand());  //impact parameter
  G4double z              = b;                            
  G4double x              = - sqrt(sqr(radius) - sqr(z));
  G4double y              = 0;

  G4ThreeVector point(x,y,z);
  return point;
}

G4ThreeVector G4NucleusModel::GetInteractionPoint(){
  return interactionPoint;
}

G4double G4NucleusModel::GetInteractionLength(){
  
   G4String particleName = incidentParticle->GetDefinition()->GetParticleName();  
   G4BertiniData* data = xSecDataBase->Instance();   
   G4double protonXSec, neutronXSec;  
   G4double incidentEnergy = incidentParticle->GetKineticEnergy();

    if(particleName == "proton"){ 
      protonXSec = data->GetCrossSection(G4ProtonProtonTotal, incidentEnergy)*1E+26; 
      neutronXSec = data->GetCrossSection(G4NeutronProtonTotal, incidentEnergy)*1E+26; 
      } 
    else{ 
        protonXSec = 1E+26; 
        neutronXSec = 0; 
    }

  G4double meanFreePath = 4*pi*myA*pow(radius0,3)/
    (3*(myZ*protonXSec + (myA - myZ)*neutronXSec)); 
  G4double lambda = - meanFreePath*log(G4UniformRand());
  return lambda;
} 

G4ThreeVector G4NucleusModel::CalculateInteractionPoint(){
  G4LorentzVector incidentMomentum = incidentParticle->Get4Momentum();
  G4double pathLength = 0.0;
  G4ThreeVector crossingPoint;
  G4ThreeVector currentPoint = incidencePoint;
  G4LorentzVector currentMomentum = incidentMomentum;
  G4VRegionModel* currentRegion = GetRegion(incidencePoint.mag());
  G4String particleName = incidentParticle->GetDefinition()->GetParticleName();
  G4ThreeVector lengthToCrossingPoint;
  G4int i;


  if(particleName == "neutron") //only protons reflect from potential wall 
    for(i=0; i<3; i++) currentPoint[i] += (currentMomentum[i]/currentMomentum.mag())*interactionLength;
    
  else{ //protons
    while(pathLength < interactionLength){
      if(!BoundaryCrossed(currentRegion, currentPoint, incidentParticle, (interactionLength-pathLength)) || currentRegion == NULL){
	for(i=0; i<3; i++) currentPoint[i] += (currentMomentum[i]/currentMomentum.mag())*(interactionLength-pathLength);
	pathLength = interactionLength;

      }
      else{
	crossingPoint = GetCrossingPoint(currentRegion, currentPoint, incidentParticle);
	currentRegion = DoBoundaryTransition(currentRegion, crossingPoint, incidentParticle);
	if(currentRegion!=NULL)
	lengthToCrossingPoint = crossingPoint - currentPoint;
	pathLength += lengthToCrossingPoint.mag();
	currentPoint = crossingPoint;
	currentMomentum = incidentParticle->Get4Momentum();
      }
    }
  }
  
   return currentPoint;
}

G4bool G4NucleusModel::BoundaryCrossed(G4VRegionModel* initialRegion,
				       G4ThreeVector initialPoint,
				       G4DynamicParticle* particle,
				       G4double distance){

  G4ThreeVector incidentMomentum = particle->GetMomentum();
  G4ThreeVector finalPoint = initialPoint;
  for(int i=0; i<3; i++) finalPoint[i] += (incidentMomentum[i]/incidentMomentum.mag())*distance;

  G4VRegionModel* finalRegion = GetRegion(finalPoint.mag());
  if(finalRegion == initialRegion) return 0;
  else return 1;
}

G4VRegionModel* G4NucleusModel::GetNextRegion(G4VRegionModel* currentRegion, G4ThreeVector initCoordinates, G4DynamicParticle* particle){
  G4ThreeVector currentMomentum = particle->GetMomentum();
  G4int i;
  G4VRegionModel* nextRegion;

  regionIterator iter;
  regionIterator begin = regions.begin();


  for(iter = regions.begin(); iter < regions.end(); iter++){
    if((*iter) == currentRegion) break;
  }

  if(currentRegion == (*begin)){
    movingIn = false;
    iter = begin + 1;
    nextRegion = *(iter);
    return nextRegion;
  } 
  
  G4double dotProduct = 0.0;

  for(i = 0; i  < 3; i++) dotProduct += initCoordinates[i]*currentMomentum[i];
  
  G4double cosAlpha = dotProduct / (currentMomentum.mag()*initCoordinates.mag());
  if(cosAlpha < -1) cosAlpha = -1;
  if(cosAlpha > 1)  cosAlpha = 1;
  G4double alpha = acos(cosAlpha);
  G4double beta = pi - alpha;
  G4double minimumDistance; //minimum distance to the center of nucleus

  if(alpha < pi/2) minimumDistance = initCoordinates.mag();
  else minimumDistance = initCoordinates.mag()*sin(beta);
  
  
  regionIterator previousRegionIter = iter-1;
  if(minimumDistance < (*previousRegionIter)->GetOuterRadius() && dotProduct < 0)
    movingIn = true; //cuts the inner region boundary 
  else movingIn = false;

   regionIterator lastRegionIter = regions.end() - 1;
 
  if(!movingIn && iter!=lastRegionIter){
    iter++;
    nextRegion = *iter;
  }
  else{
    if(iter == lastRegionIter && !movingIn){
     nextRegion = NULL;
    }
  }
  if(movingIn) {
    iter--;
    nextRegion = *iter;
  }
  return nextRegion;
}


G4VRegionModel* G4NucleusModel::DoBoundaryTransition(G4VRegionModel* currentRegion,
						     G4ThreeVector position,
						     G4DynamicParticle* particle){

  G4LorentzVector particle4Momentum = particle->Get4Momentum();
  G4double radialMomentum = 0.0;
  G4double r = position.mag();
  G4double protonMass = G4Proton::Proton()->GetPDGMass();
  G4double neutronMass = G4Neutron::Neutron()->GetPDGMass();

  G4int i;
  for(i=0; i < 3; i++) radialMomentum += position[i]*particle4Momentum[i];
  
  radialMomentum = radialMomentum/ r;

  G4VRegionModel* nextRegion = GetNextRegion(currentRegion, position, particle); 

  G4double potentialDifference;
  if(nextRegion!=NULL)  potentialDifference = nextRegion->GetProtonPotentialEnergy()- currentRegion->GetProtonPotentialEnergy(); 
  else potentialDifference = -currentRegion->GetProtonPotentialEnergy(); 
    
  G4double qv = sqr(potentialDifference) - 2.0*potentialDifference*particle4Momentum[3] + sqr(radialMomentum);
  
  for(regionIterator iter = regions.begin(); iter < regions.end(); iter++)
    if((*iter) == currentRegion) break;

  G4double p1r;
  if(qv <= 0.){ //reflection
    p1r = -radialMomentum;
    nextRegion = currentRegion;
  }

  else{ //transition to next region
    p1r = sqrt(qv);

    if (radialMomentum < 0.) p1r = -p1r;
    if(movingIn){
      iter--; nextRegion = *iter;
       }
    else
    if(nextRegion!=NULL){
      iter++; nextRegion = *iter;
    }
   }

  G4double prr = (p1r-radialMomentum)/r;
  for(i = 0; i < 3; i++) particle4Momentum[i] += position[i]*prr; //update momentum
  
  G4double totalEnergy; 
  G4double momentumMag2 = 0.0; 
  for(i=0; i < 3; i++)  momentumMag2 += sqr(particle4Momentum[i]);
  totalEnergy = sqrt(sqr(protonMass) + momentumMag2); 
  particle4Momentum[3] = totalEnergy; //fourth component
  particle->Set4Momentum(particle4Momentum);

  return nextRegion;
}

G4ThreeVector G4NucleusModel::GetCrossingPoint(G4VRegionModel* currentRegion,
					       G4ThreeVector initCoordinates,
					       G4DynamicParticle* particle){

  G4VRegionModel* nextRegion = GetNextRegion(currentRegion, initCoordinates, particle);
  regionIterator begin = regions.begin(); 

  G4ThreeVector p = particle->GetMomentum();
  G4double px = p[0], py = p[1], pz = p[2];
  G4double dotProduct = 0.0;
  G4double d; //distance to the crossing point

  for(int i = 0; i  < 3; i++) dotProduct += initCoordinates[i]*p[i];
  if(initCoordinates.mag() == 0) d = currentRegion->GetOuterRadius();
 
  else{
    G4double cosAlpha = dotProduct / (p.mag()*initCoordinates.mag()); //cos(p,r) = p dot r /(|p||r|)
    if(cosAlpha < -1) cosAlpha = -1; //possible rounding errors
    if(cosAlpha > 1)  cosAlpha = 1;
    G4double alpha = acos(cosAlpha); //the angle between p and r
    G4double beta = pi - alpha; 
    G4double cosBeta = cos(beta);
    G4double a; // third side of the triangle

    if(movingIn) a = nextRegion->GetOuterRadius();
    else a = currentRegion->GetOuterRadius();

    //the distance to crossing point is calculated from cosine formula:
    // a^2 = r^2 + d^2 -2*r*d*cos(r,d) => d = ...
 
    if(movingIn){
     	d = (2*initCoordinates.mag()*cosBeta - sqrt(4*sqr(initCoordinates.mag()*cosBeta)
						  + 4*sqr(a)-4*sqr(initCoordinates.mag())))/2;
    }
    else{
 	d = (2*initCoordinates.mag()*cosBeta + sqrt(4*sqr(initCoordinates.mag()*cosBeta)
						  + 4*sqr(a)-4*sqr(initCoordinates.mag())))/2;
    }
  }
  //let's move to the crossing point
  G4double finalX = initCoordinates[0] + (px / p.mag())*d;
  G4double finalY = initCoordinates[1] + (py / p.mag())*d;
  G4double finalZ = initCoordinates[2] + (pz / p.mag())*d;
  G4ThreeVector crossingPoint(finalX, finalY, finalZ);
  return crossingPoint;
  
}

void G4NucleusModel::DoCollision(G4DynamicParticle* particle){
  incidentParticle = particle;
  incidencePoint = GetIncidencePoint();
  interactionLength = GetInteractionLength();
  interactionPoint = CalculateInteractionPoint();
}


void G4NucleusModel::DoCollision(G4DynamicParticle* particle, G4ThreeVector coordinates){
  incidentParticle = particle;
  incidencePoint = coordinates;
  interactionLength = GetInteractionLength();
  interactionPoint = CalculateInteractionPoint();
}


G4bool G4NucleusModel::IsInside(G4ThreeVector position){
  if(position.mag() > radius) return 0;
  else return 1;
}

G4VRegionModel* G4NucleusModel::GetRegion(G4double r){
  regionIterator iter;

  for(iter =  regions.begin(); iter < regions.end(); iter++){
     if(r <= ((*iter)->GetOuterRadius())) return (*iter);
  }
  return NULL; //out of nucleus
}




G4bool G4NucleusModel::CheckPauliPrinciple(G4DynamicParticle* productParticle){
  G4double rad = interactionPoint.mag();
  G4String particleName = productParticle->GetDefinition()->GetParticleName();
  G4VRegionModel* thisRegion = GetRegion(rad);
  G4double qMax = 0;

  if(particleName == "proton")  qMax = thisRegion->GetProtonMaximumMomentum();
  if(particleName == "neutron") qMax = thisRegion->GetNeutronMaximumMomentum();
  
  G4double fermiEnergy = pow(qMax,2)/(2*productParticle->GetMass());
  G4double productEnergy = productParticle->GetKineticEnergy();

  if(productEnergy > fermiEnergy) return 1;
  else return 0;
}

G4DynamicParticle* G4NucleusModel::ReturnTargetNucleon(){ 

  //method from Methods in Subnuclear Physics
  G4DynamicParticle *targetParticle = new G4DynamicParticle();
  G4double impactRadius = interactionPoint.mag(); 
  G4VRegionModel* thisRegion = GetRegion(impactRadius); 

  G4String particleName = incidentParticle->GetDefinition()->GetParticleName();  
  G4double incidentEnergy = incidentParticle->GetKineticEnergy();

   G4BertiniData* data = xSecDataBase->Instance();   
   G4double protonXSec, neutronXSec;  

    if(particleName == "proton"){ 
      protonXSec = data->GetCrossSection(G4ProtonProtonTotal, incidentEnergy); 
      neutronXSec = data->GetCrossSection(G4NeutronProtonTotal, incidentEnergy); 
      } 

    else{ //impact particle can currently be only proton
        protonXSec = 1; 
        neutronXSec = 0; 
    }
  
    G4double ratio = myZ*protonXSec/( myZ*protonXSec + (myA - myZ)*neutronXSec); 
    G4double randomNumber = G4UniformRand();

    if(randomNumber < ratio) targetParticle = thisRegion->GenerateProton();
    if(randomNumber >= ratio) targetParticle = thisRegion->GenerateNeutron();
    
  return targetParticle; 
}

void G4NucleusModel::PrintModel(){
  cout << "A: " << myA << endl;
  cout << "Z: " << myZ << endl;
  cout << "Radius: " << radius <<" fm"<< endl;
  cout << "Atomic mass: " << atomicMass <<" MeV " << endl;
  cout << "Binding energy: " << bindingEnergy <<" MeV " << endl;
  cout << "Excitation energy: " << excitationEnergy <<" MeV " << endl;
  cout << "Kinetic energy: " << kineticEnergy <<" MeV "<< endl;
  cout << "Momentum: " << momentum <<" MeV " << endl;
   cout << "Number of regions: " << NUMBER_OF_REGIONS << endl;
  cout <<"--------------------------------------------------" << endl;
  for(regionIterator i = regions.begin(); i < regions.end(); i++){
    (*i)->PrintRegion();
    cout <<"--------------------------------------------------" << endl;
  } 
}
  
void G4NucleusModel::ApplyRecoil(trackVector track){
  G4double nucleonAvMass = 931.162;
  G4LorentzVector fourMomentum(0.0, 0.0, 0.0, 0.0);
  trackIterator iter;

  if(myA > 4){
    for(iter = track.begin(); iter < track.end(); iter++){
      fourMomentum += (*iter)->Get4Momentum();
    }

    G4ThreeVector recoilMomentum(-fourMomentum[0], -fourMomentum[1], -fourMomentum[2]);
    momentum += recoilMomentum;
  
    G4double mass = nucleonAvMass*myA;
    G4double recoilEnergy = sqrt(sqr(mass) + sqr(momentum.mag())) - mass;
    kineticEnergy = recoilEnergy; 
  }
}

/*
G4double G4NucleusModel::CalculateExcitationEnergy()
  {
    G4double incidentKineticEnergy = incidentParticle->GetKineticEnergy()/1000;
    G4double neutronMass  = G4Neutron::Neutron()->GetPDGMass() / 1000;
    G4double electronMass = 0.000511;
    G4double exnu         = 0.; 
    G4double excitationEnergyGPN   = 0.;
    G4double excitationEnergyDTA   = 0.;


    if (atomicMass/1000 > (neutronMass + 2.*electronMass))
       {
         G4int    magic = ((G4int)(myZ+0.1) == 82) ? 1 : 0;
         G4double ekin  = G4std::min(G4std::max(incidentKineticEnergy, 0.1), 4.);
         G4double cfa   = G4std::max(0.35 +((0.35 - 0.05)/2.3)*log(ekin), 0.15);
                  exnu  = 7.716*cfa*exp(-cfa);
         G4double atno  = G4std::min(atomicMass/1000, 120.);
                  cfa   = ((atno - 1.)/120.) * exp(-(atno-1.)/120.);
                  exnu  = exnu * cfa;
         G4double fpdiv = G4std::max(1.-0.25*ekin*ekin, 0.5);
         G4double gfa   = 2.*((atomicMass/1000-1.)/70.) 
                            * exp(-(atomicMass/1000-1.)/70.);

         excitationEnergyGPN = exnu * fpdiv;
         excitationEnergyDTA = exnu - excitationEnergyGPN;  
        
         G4double ran1 = 0., ran2 = 0.;
	 if (!magic)
	    { ran1 = normal();
              ran2 = normal();
            }
         excitationEnergyGPN = G4std::max(excitationEnergyGPN*(1.+ran1*gfa),0.);
         excitationEnergyDTA = G4std::max(excitationEnergyDTA*(1.+ran2*gfa),0.);
	 //cout <<"excitationEnergyGPN1: " << excitationEnergyGPN << endl;
	 //cout <<"excitationEnergyDTA1: " << excitationEnergyDTA << endl; 

         exnu = excitationEnergyGPN + excitationEnergyDTA;
	 
         if(verboseLevel > 1)
           {
             G4cout << " NuclearExcitation: " << magic << " " <<  ekin << " " << cfa
                  << " " <<  atno << " " << fpdiv << " " <<  gfa << " " << excitationEnergyGPN
                  << " " <<  excitationEnergyDTA << G4endl;
           } 
	 
         while (exnu >= incidentKineticEnergy)
	    {
              excitationEnergyGPN *= (1. - 0.5*normal());
              excitationEnergyDTA *= (1. - 0.5*normal());
	      // cout <<"excitationEnergyGPN2: " << excitationEnergyGPN << endl;
	      // cout <<"excitationEnergyDTA2: " << excitationEnergyDTA << endl; 
              exnu = excitationEnergyGPN + excitationEnergyDTA;
            }
       }             
    return exnu*1000;
}     



G4double G4NucleusModel::normal()
 {
   G4double ran = -6.0;
   for(G4int i=0; i<12; i++) ran += G4UniformRand();
   return ran;
 }

void G4NucleusModel::ApplyExcitationEnergy(){
  excitationEnergy += CalculateExcitationEnergy();
}

 */















