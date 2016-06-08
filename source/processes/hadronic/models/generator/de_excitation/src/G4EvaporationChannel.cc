// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
// some corrections V.Krylov
//

#include "G4EvaporationChannel.hh"

G4EvaporationChannel::G4EvaporationChannel(const G4int theGamma,
					   const G4int theA,
					   const G4int theZ,
					   RWTValVector<G4double> * theExcitationEnergies,
					   RWTValVector<G4int> * theExcitationSpins):
  Gamma(theGamma),
  A(theA),
  Z(theZ),
  ExcitationEnergies(theExcitationEnergies),
  ExcitationSpins(theExcitationSpins),
  AResidual(0),
  ZResidual(0),
  CoulombBarrier(0.0),
  BindingEnergy(0.0),
  MaximalKineticEnergy(-1000.0),
  EmissionProbability(0.0)
{ 
  theEvaporationProbabilityPtr = new G4EvaporationProbability(this);
  MyOwnEvaporationProbability = true;

  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
  MyOwnLevelDensity = true;
}

G4EvaporationChannel::~G4EvaporationChannel()
{

  if (MyOwnEvaporationProbability) delete theEvaporationProbabilityPtr;

  if (MyOwnLevelDensity) delete theLevelDensityPtr;
}




G4EvaporationChannel::G4EvaporationChannel(const G4EvaporationChannel & right)
{
  G4Exception("G4EvaporationChannel::copy_costructor meant to not be accessable");
}

const G4EvaporationChannel & G4EvaporationChannel::operator=(const G4EvaporationChannel & right)
{
  G4Exception("G4EvaporationChannel::operator= meant to not be accessable");
  return *this;
}

G4bool G4EvaporationChannel::operator==(const G4EvaporationChannel & right) const 
{
  return (this == (G4EvaporationChannel *) &right);
    //  return false;
}

G4bool G4EvaporationChannel::operator!=(const G4EvaporationChannel & right) const 
{
    return (this != (G4EvaporationChannel *) &right);
    //  return true;
}



void G4EvaporationChannel::Initialize(const G4Fragment & fragment)
{

  G4int anA = fragment.GetA();
  G4int aZ = fragment.GetZ();
  G4double ExEnergy = fragment.GetExcitationEnergy();

  AResidual = anA - A;
  ZResidual = aZ - Z;

  // We only take into account channels which are physically allowed
  if (AResidual <= 0 || ZResidual <= 0 || AResidual < ZResidual ||
      (AResidual == ZResidual && AResidual > 1)) {
    LevelDensityParameter = 0.0;
    CoulombBarrier = 0.0;
    BindingEnergy = 0.0;
    MaximalKineticEnergy = -1000.0*MeV;
    EmissionProbability = 0.0;
  } else {
    // Get Level Density
    LevelDensityParameter = theLevelDensityPtr->LevelDensityParameter(anA,aZ,ExEnergy);

    // Coulomb Barrier calculation
    CoulombBarrier = CalcCoulombBarrier(AResidual,ZResidual)*MeV;
	
    
    // Binding Enegy (for separate fragment from nucleus)
    BindingEnergy = CalcBindingEnergy(anA,aZ)*MeV;

    // Maximal Kinetic Energy
    MaximalKineticEnergy = CalcMaximalKineticEnergy(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ,anA)+ExEnergy)*MeV;

    // Emission probability
    if (MaximalKineticEnergy <= 0.0) EmissionProbability = 0.0;
    else { 
      // Total emission probability for this channel
      EmissionProbability = theEvaporationProbabilityPtr->EmissionProbability(fragment,0.0);
      // Next is a loop over excited states for this channel summing probabilities
      G4double SavedGamma = Gamma;
      G4double SavedMaximalKineticEnergy = MaximalKineticEnergy;
      for (G4int i = 0; i < ExcitationEnergies->length(); i++) {
	if (ExcitationSpins->operator()(i) < 0.1) continue;
	Gamma = ExcitationSpins->operator()(i)*A;
	// substract excitation energies
	MaximalKineticEnergy -= ExcitationEnergies->operator()(i)/MeV;
	// update probability
	G4double tmp = theEvaporationProbabilityPtr->EmissionProbability(fragment,0.0);
	EmissionProbability += tmp;
      }
      // restore Gamma and MaximalKineticEnergy
      MaximalKineticEnergy = SavedMaximalKineticEnergy;
      Gamma = SavedGamma;
    }
  }

  return;

}


G4FragmentVector * G4EvaporationChannel::BreakUp(const G4Fragment & theNucleus)
{
  // calculate kinetic energy of evaporated fragment
  G4double EvaporatedKineticEnergy = CalcKineticEnergy();  // MeV
  G4double EvaporatedMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A)/MeV; // MeV
  G4double EvaporatedEnergy = EvaporatedKineticEnergy + EvaporatedMass;


  G4ThreeVector momentum( IsotropicVector( sqrt( EvaporatedEnergy*EvaporatedEnergy - 
						 EvaporatedMass*EvaporatedMass )
					   ) );

  G4LorentzVector EvaporatedMomentum( momentum, EvaporatedEnergy );
  EvaporatedMomentum.boost( theNucleus.GetMomentum().boostVector() );
  
  // to avoid rounding errors in Lorentz boost which then produce 
  // evaporated fragments with excitation energies ~10^-10 eV
  EvaporatedMomentum.setE(sqrt(EvaporatedMomentum.vect().mag2()+EvaporatedMass*EvaporatedMass));
  
  
  

  G4Fragment * EvaporatedFragment = new G4Fragment( A, Z, EvaporatedMomentum );
  if ( !EvaporatedFragment )
    G4Exception( "G4EvaporationChannel::BreakUp: Can't create G4Fragment! ");

  G4LorentzVector FragmentMomentum( theNucleus.GetMomentum() );
  FragmentMomentum.boost( -theNucleus.GetMomentum().boostVector() );
  

  G4LorentzVector ResidualMomentum( -momentum, FragmentMomentum.e() - EvaporatedEnergy );
  ResidualMomentum.boost( theNucleus.GetMomentum().boostVector() );

  
  G4Fragment * ResidualFragment = new G4Fragment( AResidual, ZResidual, ResidualMomentum );
  if ( !ResidualFragment )
		G4Exception( "G4EvaporationChannel::BreakUp: Can't create G4Fragment! ");


  G4FragmentVector * theResult = new G4FragmentVector;
  if ( !theResult )
    G4Exception( "G4EvaporationChannel::BreakUp: Can't create G4FragmentVector! ");

  theResult->insert(EvaporatedFragment);
  theResult->insert(ResidualFragment);

  return theResult; 
} 


G4double G4EvaporationChannel::CalcCoulombBarrier(const G4int ARes, const G4int ZRes)
  // Calculation of Coulomb potential energy (barrier) in MeV for outgoing fragment
{
  G4double Barrier = 0.0;
  if (Z == 0 && A == 1) return 0.0; // for neutron
  else {
    G4int nZZRes = Z * ZRes;
    G4double r0 = 2.173*(1.0+0.006103 * nZZRes)/(1.0+0.009443 * nZZRes);
    Barrier = 1.44/r0 * nZZRes / (pow( A,1./3. ) + pow( ARes,1./3. ));
  }
  return Barrier;
}



G4double G4EvaporationChannel::CalcBindingEnergy(const G4int anA, const G4int aZ)
  // Calculate Binding Energy for separate fragment from nucleus
{
  // Mass Excess for residual nucleus
  G4double ResNucMassExcess = G4NucleiProperties::GetMassExcess(AResidual,ZResidual)/MeV;
  // Mass Excess for fragment
  G4double FragmentMassExcess = G4NucleiProperties::GetMassExcess(A,Z)/MeV;
  // Mass Excess for Nucleus
  G4double NucleusMassExcess = G4NucleiProperties::GetMassExcess(anA,aZ)/MeV;

  return ResNucMassExcess + FragmentMassExcess - NucleusMassExcess;
}


G4double G4EvaporationChannel::CalcMaximalKineticEnergy(const G4double NucleusTotalE)
  // Calculate maximal kinetic energy that can be carried by fragment (in MeV)
{
//   // Odd-Even term correction (for maximal kinetic energy)
//   G4double odd = 0.0;
//   if (A < 65) {
//     G4int NCorr = A - Z;
//     NCorr = 2*(NCorr/2) - NCorr;
//     G4int ZCorr = 2*(Z/2) - Z;
//     odd = 11.0*(2+NCorr+ZCorr)/sqrt(G4double(A));
//   }

//   if (A <= 55) return U/MeV - 
// 		 (BindingEnergy + CoulombBarrier)/MeV - odd;
//   else if (A > 55 && A < 65) return U/MeV - 
// 			       (BindingEnergy + CoulombBarrier)/MeV -
// 			       odd * (1.0 - (A-55)/10.);

//   else return U/MeV - (BindingEnergy + CoulombBarrier);

  G4double ResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( ZResidual, AResidual )/MeV;
  G4double EvaporatedMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( Z, A )/MeV;

  return ( (NucleusTotalE/MeV)*(NucleusTotalE/MeV) + 
	   EvaporatedMass*EvaporatedMass - ResidualMass*ResidualMass)/
    (2.0*NucleusTotalE  ) -
    EvaporatedMass - CoulombBarrier/MeV;
}




G4double G4EvaporationChannel::CalcKineticEnergy(void)
  // Samples fragment kinetic energy (in MeV).
  // It uses Dostrovsky's approximation for the inverse reaction cross
  // in the probability for fragment emisson
{
  if (MaximalKineticEnergy < 0.0) 
    G4Exception("G4EvaporationChannel::CalcKineticEnergy: maximal kinetic energy is less than 0");
  
//   G4double Rb = 4.0*LevelDensityParameter/(1./MeV)*AResidual*MaximalKineticEnergy/MeV;
//   G4double RbSqrt = sqrt(Rb);
//   G4double PEX1 = 0.0;
//   if (RbSqrt < 160.0) PEX1 = exp(-RbSqrt);
//   G4double Rk = 0.0;
//   G4double FRk = 0.0;
//   do {
//     G4double RandNumber = G4UniformRand();
//     Rk = 1.0 + (1./RbSqrt)*log(RandNumber + (1.0-RandNumber)*PEX1);
//     G4double Q1 = 1.0;
//     G4double Q2 = 1.0;
//     if (Z == 0) { // for emitted neutron
//       G4double Beta = (2.12/pow(AResidual,2./3.) - 0.05)/
// 	(0.76 + 2.2/pow(AResidual,1./3.));
//       Q1 = 1.0 + Beta/(MaximalKineticEnergy/MeV);
//       Q2 = Q1*sqrt(Q1);
//     } 
    
//     FRk = (3.0*sqrt(3.0)/2.0)/Q2 * Rk * (Q1 - Rk*Rk);
    
//   } while (FRk < G4UniformRand());


  G4double Rb = 4.0*LevelDensityParameter/(1./MeV)*AResidual*(MaximalKineticEnergy)/MeV;
  G4double RbSqrt = sqrt(Rb);
  G4double PEX1 = 0.0;
  if (RbSqrt < 160.0) PEX1 = exp(-RbSqrt);
  G4double Rk = 0.0;
  G4double FRk = 0.0;
  do {
    G4double RandNumber = G4UniformRand();
    Rk = 1.0 + (1./RbSqrt)*log(RandNumber + (1.0-RandNumber)*PEX1);
    G4double Q1 = 1.0;
    G4double Q2 = 1.0;
    if (Z == 0) { // for emitted neutron
      G4double Beta = (2.12/pow(AResidual,2./3.) - 0.05)/
	(0.76 + 2.2/pow(AResidual,1./3.));
      Q1 = 1.0 + Beta/(MaximalKineticEnergy/MeV);
      Q2 = Q1*sqrt(Q1);
    } 
    
    FRk = (3.0*sqrt(3.0)/2.0)/Q2 * Rk * (Q1 - Rk*Rk);
    
  } while (FRk < G4UniformRand());



  G4double result =  (MaximalKineticEnergy)/MeV * (1.0-Rk*Rk) + CoulombBarrier/MeV;

  return result;
} 
 

G4ThreeVector G4EvaporationChannel::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*cos(Phi)*SinTheta,
                       Magnitude*sin(Phi)*SinTheta,
                       Magnitude*CosTheta);
  return Vector;
}



