#include "../include/G4NucleusModel.hh"

G4NucleusModel::G4NucleusModel(const G4int A, const G4int Z){
  massNumber = A;
  protonNumber = Z;
  //bindingEnergy = 7*MeV;
  excitationEnergy = 0.0;
  bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z,A);
  shellStructure = new G4RegionModel(NUMBER_OF_REGIONS, A, Z);   
  radius = radius0*pow(massNumber, 1./3.);
}

G4NucleusModel::~G4NucleusModel(){
  delete shellStructure;
}

void G4NucleusModel::SetParameters(const G4int A, const G4int Z){
  massNumber = A;
  protonNumber = Z;
  bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z,A);
}

void G4NucleusModel::AddExcitationEnergy(G4double additionalEnergy){
  excitationEnergy += additionalEnergy;
}

G4double G4NucleusModel::GetMassNumber(){
  return massNumber;
}

G4double G4NucleusModel::GetProtonNumber(){
  return protonNumber;
}

G4double G4NucleusModel::GetNuclearRadius(){
  return radius;
}

G4double G4NucleusModel::GetExcitationEnergy(){
  return excitationEnergy;
}

G4double G4NucleusModel::GetNuclearMass(){
   return protonNumber*G4Proton::Proton()->GetPDGMass() + 
          (massNumber-protonNumber)*G4Neutron::Neutron()->GetPDGMass() -
          bindingEnergy;
}

G4double G4NucleusModel::GetAtomicMass(){


    const G4double electron_mass = G4Electron::Electron()->GetPDGMass()/MeV;
    const G4double proton_mass = G4Proton::Proton()->GetPDGMass()/MeV;
    const G4double neutron_mass = G4Neutron::Neutron()->GetPDGMass()/MeV;
    //
    // Weitzsaecker's Mass formula
    //
    G4int myZ = G4int(protonNumber + 0.5);
    G4int myA = G4int(massNumber + 0.5);
      
    G4double mass =
      (massNumber-protonNumber)*neutron_mass + protonNumber*proton_mass + protonNumber*electron_mass
      - 15.67*massNumber                                          // nuclear volume
      + 17.23*pow(massNumber,2./3.)                               // surface energy
      + 93.15*pow(massNumber/2.-protonNumber,2.)/massNumber       // asymmetry
      + 0.6984523*pow(protonNumber,2.)*pow(massNumber,-1./3.);    // coulomb
    G4int ipp = (myA - myZ)%2;                                    // pairing
    G4int izz = myZ%2;
    if( ipp == izz )mass += (ipp+izz-1) * 12.0 * pow(massNumber,-0.5);
    return mass*MeV;
}

G4double G4NucleusModel::GetNuclearDensity(G4ThreeVector position){
  G4double r=position.mag();
  return (shellStructure->GetDensity(r));
}

void G4NucleusModel::DoCollision(G4ThreeVector initialMomentum, G4ParticleType particle){
  incidentMomentum = initialMomentum;
  incidentParticle = particle;
  //G4double pmass = particle.mass;
  //incidentEnergy = pow(initialMomentum.mag(), 2)/(2*pmass);
}

G4ThreeVector G4NucleusModel::GetIncidencePoint(){

  G4double b              = radius * sqrt(G4UniformRand());  //impact parameter
  cout << "Impact parameter: " << b << G4endl;
  G4double z              = b;                            
  G4double x              = - sqrt(sqr(radius) - sqr(z));
  G4double y              = 0;

  G4ThreeVector point(x,y,z);
  //incidencePoint = point;
  return point;
}

G4double G4NucleusModel::GetInterActionLength(){
  /*  
    G4HETCData data;
    switch(incidentParticle){
    case G4Proton:
      protonXSec = data.GetCrossSection(ProtonProtonTotal, incidentEnergy);
      neutronXSec = data.GetCrossSection(NeutronProtonTotal, incidentEnergy);
      break;
    } 
    meanFreePath = 4*pi*massNumber*pow(radius0,3)/
    (3*(protonNumber*protonXSec + (massNumber - protonNumber)*neutronXSec));
  */   

  G4double meanFreePath = 0.5E-15;
  
  G4double lambda = - meanFreePath*log(G4UniformRand());
  cout << "Interaction length: " << lambda << G4endl;
  return lambda;
} 

G4ThreeVector G4NucleusModel::GetInterActionPoint(){
  incidencePoint = GetIncidencePoint();
  interActionPoint = incidencePoint;
  interActionLength = GetInterActionLength();

  interActionPoint[0] += interActionLength;
  return interActionPoint;
}

G4bool G4NucleusModel::IsInside(G4ThreeVector position){
  if(position.mag() > radius) return 0;
  else return 1;
}

G4bool G4NucleusModel::PauliPrinciple(G4double productEnergy, G4double productMass){
  G4double qMax = shellStructure->GetMaximumNucleonMomentum(0, 0);
  //G4double qMax = 5; 
  if(productEnergy > pow(qMax,2)/(2*productMass)) return 0;
  else return 1;
}

G4ParticleType G4NucleusModel::ReturnTargetNucleon(){  //does not work!

  /*
    G4HETCData data;
    switch(incidentParticle){
    case G4Proton:
      protonXSec = data.GetCrossSection(ProtonProtonTotal, incidentEnergy);
      neutronXSec = data.GetCrossSection(NeutronProtonTotal, incidentEnergy);
      break;
    default:
       protonXSec = 1;
       neutronXSec = 0;
    }

    G4double ratio = protonNumber*protonXSec/( protonNumber*protonXSec + (massNumber - protonNumber)*neutronXSec);
    G4double randomNumber = G4UniformRand();
    if(randomNumber < ratio)
      return G4Proton;
    else
      return G4Neutron;
  */  

  targetParticle = G4Proton;
  return targetParticle;
}

G4ThreeVector G4NucleusModel::GetTargetMomentum(){
  //choose direction and magnitude
  G4double cosMu  = 1 - 2 * G4UniformRand();
  G4double phi = 2 * pi * G4UniformRand();
  G4double mu = acos(cosMu);

  G4double qMax = shellStructure->GetMaximumNucleonMomentum(0, 0);
  cout << "qMax: " << qMax << G4endl;
  //G4double qMax = 5;
  G4double q = qMax * sqrt(G4UniformRand());
  G4ThreeVector p(q*cos(mu)*sin(phi), q*sin(mu)*sin(phi), q*cos(phi));

  targetMomentum = p;
  return targetMomentum;
}  

/*
G4NuclearFermiDensity::G4NuclearFermiDensity(G4double anA, G4double aZ) 
  :  a(0.545 * fermi) 
//        const G4double r0=1.14*fermi;
	const G4double r0=1.16 * ( 1. - 1.16 * pow(anA, -2./3.)) * fermi;
	theA = G4int(anA);
	theZ = G4int(aZ);
	theR= r0 * pow(theA, 1./3. );
	Setrho0(3./ (4. * pi * pow(r0,3.) * theA * ( 1. + sqr(a/theR)*pi2 )));
class G4FermiMomentum :
G4double cbrt(G4double x) { return pow(x,1./3.); }
inline G4double GetFermiMomentum(G4double density)
    {
	return constofpmax * cbrt(density * theA);
   }

constofpmax(hbarc*cbrt(3.*pi2))

inline G4ThreeVector GetMomentum(G4double density)
    { 
	G4ThreeVector p;
	
	do {
	    p=G4ThreeVector(2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.,
	    		    2.*G4UniformRand()-1.);
	    } while ( p.mag() > 1. );
	return p*GetFermiMomentum(density); 
     }

G4double
G4HEInelastic::NuclearExcitation(G4double  incidentKineticEnergy,
                                 G4double  atomicWeight,
                                 G4double  atomicNumber,
                                 G4double& excitationEnergyGPN,
                                 G4double& excitationEnergyDTA)
  {
    G4double neutronMass  = Neutron.getMass();
    G4double electronMass = 0.000511;
    G4double exnu         = 0.; 
    excitationEnergyGPN   = 0.;
    excitationEnergyDTA   = 0.;

    if (atomicWeight > (neutronMass + 2.*electronMass))
       {
         G4int    magic = ((G4int)(atomicNumber+0.1) == 82) ? 1 : 0;
         G4double ekin  = Amin(Amax(incidentKineticEnergy, 0.1), 4.);
         G4double cfa   = Amax(0.35 +((0.35 - 0.05)/2.3)*log(ekin), 0.15);
                  exnu  = 7.716*cfa*exp(-cfa);
         G4double atno  = Amin(atomicWeight, 120.);
                  cfa   = ((atno - 1.)/120.) * exp(-(atno-1.)/120.);
                  exnu  = exnu * cfa;
         G4double fpdiv = Amax(1.-0.25*ekin*ekin, 0.5);
         G4double gfa   = 2.*((atomicWeight-1.)/70.) 
                            * exp(-(atomicWeight-1.)/70.);

         excitationEnergyGPN = exnu * fpdiv;
         excitationEnergyDTA = exnu - excitationEnergyGPN;  
        
         G4double ran1 = 0., ran2 = 0.;
	 if (!magic)
	    { ran1 = normal();
              ran2 = normal();
            }
         excitationEnergyGPN = Amax(excitationEnergyGPN*(1.+ran1*gfa),0.);
         excitationEnergyDTA = Amax(excitationEnergyDTA*(1.+ran2*gfa),0.);
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
              exnu = excitationEnergyGPN + excitationEnergyDTA;
            }
       }             
    return exnu;

    const G4double electron_mass = G4Electron::Electron()->GetPDGMass()/MeV;
    const G4double proton_mass = G4Proton::Proton()->GetPDGMass()/MeV;
    const G4double neutron_mass = G4Neutron::Neutron()->GetPD
    //
    // Weitzsaecker's Mass formula
    //
    G4double mass =
      (A-Z)*neutron_mass + Z*proton_mass + Z*electron_mass
      - 15.67*A                                          // nuclear volume
      + 17.23*pow(A,2./3.)                               // surface energy
      + 93.15*pow(A/2.-Z,2.)/A                           // asymmetry
      + 0.6984523*pow(Z,2.)*pow(A,-1./3.);               // coulomb
    G4int ipp = (myA - myZ)%2;            // pairing
    G4int izz = myZ%2;
    if( ipp == izz )mass += (ipp+izz-1) * 12.0 * pow(A,-0.5);
    return mass*MeV;
  }

 */















