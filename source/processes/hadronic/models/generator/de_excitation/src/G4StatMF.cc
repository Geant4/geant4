
#include "G4StatMF.hh"



// Default constructor
G4StatMF::G4StatMF() : theSim(0)
  //:theMicrocanonicalSim(0),theMacrocanonicalSim(0)
{

}


// Destructor
G4StatMF::~G4StatMF() 
{
  //  if (theMicrocanonicalSim != 0) delete theMicrocanonicalSim;

  //  if (theMacrocanonicalSim != 0) delete theMacrocanonicalSim;

}


// Operators

G4StatMF & G4StatMF::operator=(const G4StatMF & right)
{
  G4Exception("G4StatMF::operator= meant to not be accessable");
  return *this;
}


G4bool G4StatMF::operator==(const G4StatMF & right)
{
  return false;
}
 

G4bool G4StatMF::operator!=(const G4StatMF & right)
{
  return true;
}





G4FragmentVector * G4StatMF::BreakItUp(const G4Fragment &theFragment)
{
  G4StatMFMicrocanonical * theMicrocanonicalSim = 0;
  G4StatMFMacrocanonical * theMacrocanonicalSim = 0;


  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 and M_0 = 3.3 for A <= 110
  G4double MaxAverageMultiplicity = 2.6;
  if (theFragment.GetA() <= 110) MaxAverageMultiplicity = 3.3;

  if (theFragment.GetExcitationEnergy()/MeV <= 0.0) return 0;

  //-------------------------------------------------------
  // The first part of simulation procedure
  // Direct simulation part (Microcanonical simulation)
  //-------------------------------------------------------
  theMicrocanonicalSim = new G4StatMFMicrocanonical(theFragment);

  G4int Iterations = 0;
  G4double Temperature = 0.0;
  G4double EnergyCoulomb = 0.0;
  
  G4bool FirstTime = true;

  do {
    G4bool StrangeFragment = false;
    do {
      
      //      if (theMicrocanonicalSim->GetMeanMultiplicity() <= MaxAverageMultiplicity) {
      G4double theMeanMult = theMicrocanonicalSim->GetMeanMultiplicity();
      if (theMeanMult <= MaxAverageMultiplicity) {

	// choose fragments atomic numbers and charges from direct simulation
	theMicrocanonicalSim->ChooseAandZ(theFragment);
	theSim = theMicrocanonicalSim;
      } else {
	//-------------------------------------------------
	// Non direct part (Macrocanonical Simulation)
	//-------------------------------------------------
	if (FirstTime) {
	  theMacrocanonicalSim = new G4StatMFMacrocanonical(theFragment);
	  theSim = theMacrocanonicalSim;
	  FirstTime = false;
	}
	// Select calculated fragment total multiplicity, 
	// fragment atomic numbers and fragment charges.
	theMacrocanonicalSim->ChooseAandZ(theFragment);
      }
      
      StrangeFragment = false;
      for (G4int i = 0; i < theSim->GetMeanMultiplicity(); i++) {
	G4int A = theSim->GetFragmentA(i);
	G4int Z = theSim->GetFragmentZ(i);
	if (A > 1 && (Z >= A || Z <= 0)) {
	  StrangeFragment = true;
	  break;
	}
      }

    } while (StrangeFragment);

    // I don't know what to do with this
    // ????????????????????????????????????
    //    if (theSim->GetMultiplicity() <= 1) Evaporation();

    //--------------------------------------
    // Second part of simulation procedure.
    //--------------------------------------

    // Find temperature of breaking channel.
    Temperature = theSim->GetMeanTemperature();

    if (FindTemperatureOfBreakingChannel(theFragment,theSim->GetMultiplicity(),
					 Temperature,EnergyCoulomb)) break;

  } while (Iterations++ < 10);


  // Separate neutrons from charged fragments
  theSim->SortFragments();


  G4ThreeVector * FragmentsMomenta = new G4ThreeVector[theSim->GetMultiplicity()];

  CoulombImpulse(theFragment,theSim->GetNumOfCharged(),theSim->GetMultiplicity(),
		 theSim->GetMeanTemperature(),EnergyCoulomb,FragmentsMomenta);


  if (theSim->GetNumOfCharged() < theSim->GetMultiplicity()) 
    CalculateFragmentsMomentum(theSim->GetNumOfNeutrons(),theSim->GetNumOfCharged(),Temperature,
			     1.5*Temperature*theSim->GetNumOfNeutrons(),FragmentsMomenta);


  // Perform Lorentz boost
  G4LorentzVector * FourMomenta = new G4LorentzVector[theSim->GetMultiplicity()];
  G4int i;
  for (i = 0; i < theSim->GetMultiplicity(); i ++) {
    FourMomenta[i].setVect(FragmentsMomenta[i]);
    FourMomenta[i].setE(sqrt(FragmentsMomenta[i].mag2()+
			     G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmentA(i)*
			     G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmentA(i)));
    FourMomenta[i].boost(theFragment.GetMomentum().boostVector());
  }

  // Calculate Fragments excitation energies
  G4double * UFragments = new G4double[theSim->GetMultiplicity()];
  for (i = 0; i < theSim->GetMultiplicity(); i ++) {
    if (G4int(theSim->GetFragmentA(i)) < 3) UFragments[i] = 0.0;
    else UFragments[i] = CalculateFragmentExcitationEnergy(i,Temperature);
  }


  // insert fragments in vector
  G4FragmentVector * theResult = new G4FragmentVector;
  for (i = 0; i < theSim->GetMultiplicity(); i++) {
    G4Fragment * afragment = new G4Fragment(theSim->GetFragmentA(i),theSim->GetFragmentZ(i),
					    FourMomenta[i]);
    afragment->SetExcitationEnergy(UFragments[i]*MeV);
    theResult->insert(afragment);
  }

  delete theMicrocanonicalSim;
  if (theMacrocanonicalSim != 0) delete theMacrocanonicalSim;
  delete [] FragmentsMomenta;
  delete [] FourMomenta;
  delete [] UFragments;


  return theResult;
}


G4bool G4StatMF::FindTemperatureOfBreakingChannel(const G4Fragment & theFragment,
						  const G4double & Multiplicity,
						  G4double & Temperature,
						  G4double & EnergyCol)
  // This finds temperature of breaking channel.
{
  G4double A = theFragment.GetA();
  G4double Z = theFragment.GetZ();
  G4double U = theFragment.GetExcitationEnergy()/MeV;

  G4double A13 = pow(A,1./3.);
  G4double PkP13 = pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.); 



  const G4double HT = 0.5;
  G4double H = 0.0;
  G4int counter = 0;
  G4int id = 0;
  G4int K1 = 0, K2 = 0;

  G4double T = max(Temperature,0.0012);

  G4double MassExcess0 = G4NucleiProperties::GetMassExcess(A,Z);

  do {
    G4double ExchangeEnergy = 0.0;
    G4double EB = -MassExcess0;
    G4double ECP = 0.0;
    G4double ETP = 1.5*T*Multiplicity;
    G4int i ;
    for (i = 0; i < Multiplicity; i++) {
      EB += G4NucleiProperties::GetMassExcess(theSim->GetFragmentA(i),theSim->GetFragmentZ(i));
      G4double PA13 = pow(theSim->GetFragmentA(i),1.0/3.0);
      ECP += (3./5.)*(1.66/G4StatMFParameters::Getr0())*
	(theSim->GetFragmentZ(i)*theSim->GetFragmentZ(i)/PA13);

      if (G4int(theSim->GetFragmentA(i)) <= 3) continue;

      G4double ESUF = theSim->GetFragmentA(i)*T*T*2.5*G4StatMFParameters::GetBeta0()/
	(G4StatMFParameters::GetCriticalTemp()*G4StatMFParameters::GetCriticalTemp()*PA13);

      if (theSim->DBetaDT(T) == 0.0) ESUF = 0.0;
      ExchangeEnergy += theSim->GetFragmentA(i)*T*T/theSim->GetFragmentInvLevelDensity(i) + ESUF;
      if (G4int(theSim->GetFragmentA(i)) == 4) ExchangeEnergy -= ESUF;
    }
    
    EnergyCol = (3./5.)*(1.44/G4StatMFParameters::Getr0())*Z*Z*PkP13/A13 - ECP*PkP13;

    G4double TotalEnergy = ETP + ExchangeEnergy + EB + EnergyCol;

    G4double D = (U - TotalEnergy)/U;
    if (abs(D) < 0.003) {
      Temperature = T;
      return false;
    }
    counter++;

    if (D < 0.0) H = -HT;
    else H = HT;

    if (D <= 0.0) {
      for (;;) {
	K1 = 1;
	if (K1 ==1 && K2 == 1) id++;
	if (id > 30) return true;
	T += H/pow(2.0,id);
	if ( T >= 0.001) break;
	K2 = 1;
	H = HT;
      }
    } else {
      for (;;) {
	K2 = 1;
	if (K1 ==1 && K2 == 1) id++;
	if (id > 30) return true;
	T += H/pow(2.0,id);
	if ( T >= 0.001) break;
	K1 = 1;
	H = HT;
      }
    }
  } while (counter <= 120);
  return true;
}




void G4StatMF::CoulombImpulse(const G4Fragment & theFragment,
			      const G4int & NumberOfChargedFragments,
			      const G4int & Multiplicity,
			      const G4double & Temperature,
			      const G4double & CoulombEnergy,
			      G4ThreeVector * MomentumOfFragments)
  // Calculate asymptotic fragments momenta (after breakup fragments
  // will fly away under Coulomb field)
{
  G4ThreeVector * Position = new G4ThreeVector[Multiplicity];
  Place(theFragment,Multiplicity,Position);

  G4double TotalKineticEnergyOfFragments = (3./2.)*NumberOfChargedFragments*Temperature;
  G4double pedo = NumberOfChargedFragments;
  G4ThreeVector * Velocities = new G4ThreeVector[NumberOfChargedFragments];
  CalculateFragmentsMomentum(NumberOfChargedFragments,0,Temperature,
			     TotalKineticEnergyOfFragments,Velocities);
  G4int i;
  for (i = 0; i < NumberOfChargedFragments; i++) 
    //    Velocities[i].setMag(Velocities[i].mag()/(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmenA(i)/MeV));
    Velocities[i] *= (G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmentA(i)/MeV);
  
  // Solve equations of motion for fragments
  SolveEqOfMotion(Position,Velocities,MomentumOfFragments,NumberOfChargedFragments,
		  CoulombEnergy,TotalKineticEnergyOfFragments);

  for (i = 0; i < NumberOfChargedFragments; i++)
    MomentumOfFragments[i] = Velocities[i]*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmentA(i)/MeV);

  delete [] Velocities;
  delete [] Position;
  
  return;
}


void G4StatMF::Place(const G4Fragment & theFragment,
		     const G4int & Multiplicity,
		     G4ThreeVector * Position)
  // Randomly samples fragments positions inside prolongated ellipsoid
{
  const G4double RN = 1.17;
  const G4double RSys = 2.0*RN*pow(theFragment.GetA(),1./3.);
  G4int i;

  for (;;) {
    i = 1;
    // This gives the position at the breakup instant
    G4double R = (RSys - RN*pow(theSim->GetFragmentA(0),1./3.))*pow(G4UniformRand(),1./3.);
    Position[0] = IsotropicVector(R);

    G4bool again = false;
    for (;;) {
      G4int k = i;
      G4int kk = 0;
      i++;
      G4double RR = 0.0;
      G4double RMin = 0.0;
      do {
	kk++;
	if (again = (kk > 1000)) break;
	R = (RSys - RN*pow(theSim->GetFragmentA(0),1./3.))*pow(G4UniformRand(),1./3.);
	Position[i] = IsotropicVector(R);

	for (G4int j = 0; j < k; j++) {
	  G4ThreeVector tmp = Position[i] - Position[j];
	  RR = tmp.mag2();
	  RMin = RN*(pow(theSim->GetFragmentA(i),1./3.)+pow(theSim->GetFragmentA(j),1./3.));
	  if (RR < RMin*RMin) break;
	}

      } while (RR < RMin*RMin);

      if (again) {
	again = false;
	break;
      }
      if (i == Multiplicity) return;
    }
  }
}



void G4StatMF::SolveEqOfMotion(G4ThreeVector * InitialPos,
			       G4ThreeVector * InitialVel,
			       G4ThreeVector * FinalVel,
			       const G4int & Multiplicity,
			       const G4double & CoulombEnergy,
			       const G4double & KineticEnergy)
  // This method will find a solution of Newton's equation of motion
  // for fragments in the self-consistent time-dependent Coulomb field
{
  if (CoulombEnergy <= 0.0) return;

  G4int Iterations = 0;
  G4double TN = 0.0;
  G4double TS = 0.0;
  G4double DT = 2.0;

  G4ThreeVector * A = new G4ThreeVector[Multiplicity];
  G4ThreeVector * Pos = new G4ThreeVector[Multiplicity];

  G4int i;
  for (i = 0; i < Multiplicity; i++) {
    FinalVel[i] = InitialVel[i];
    Pos[i] = InitialPos[i];
  }

  G4ThreeVector ** Force = new G4ThreeVector*[Multiplicity];
  for (i = 0; i < Multiplicity; i++) 
    Force[i] = new G4ThreeVector[Multiplicity];

  G4ThreeVector * ForceS = new G4ThreeVector[Multiplicity];

  G4ThreeVector * Accel = new G4ThreeVector[Multiplicity];

  G4ThreeVector * SavedVelo = new G4ThreeVector[Multiplicity];
  
  do {

    G4int i;
    G4ThreeVector distance;
    for (i = 0; i < Multiplicity; i++) {
      for (G4int j = 0; j < Multiplicity; j++) {
	if (i != j) {
	  distance = InitialPos[i] - InitialPos[j];
	  Force[i][j] = (1.44*(theSim->GetFragmentA(i)*theSim->GetFragmentA(j))/
			 (distance.mag2()*distance.mag()))*distance;
	}
      }
    }


    for ( i = 0; i < Multiplicity; i++) {
      for (G4int j = 0; j < Multiplicity; j++) {
	if (i != j) ForceS[i] += Force[i][j];
      }
    }

    for ( i = 0; i < Multiplicity; i++) {
      Accel[i] = ForceS[i];
      Accel[i] *= 1./(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)*theSim->GetFragmentA(i)/MeV);
    }

    TN = TS + DT;

    for ( i = 0; i < Multiplicity; i++) {
      SavedVelo[i] = FinalVel[i];
      FinalVel[i] = Accel[i]*(TN-TS);
      Pos[i] += (SavedVelo[i]+FinalVel[i])*(TN-TS)*0.5;
    }

    if (Iterations >= 50 && Iterations < 75) DT = 4.;
    else if (Iterations >= 75) DT = 10.;

    TS = TN;

  } while (Iterations++ < 100);
  
  // Summed fragment kinetic energy
  G4double SummedKineticEnergy = 0.0;
  for ( i = 0; i < Multiplicity; i++) SummedKineticEnergy += 
					     theSim->GetFragmentA(i)*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*
					     0.5*FinalVel[i].mag2();

  // Scaling of fragment velocities
  G4double Eta = ( CoulombEnergy + KineticEnergy ) / SummedKineticEnergy;
  for ( i = 0; i < Multiplicity; i++) FinalVel[i] *= Eta;


  // Garbage collection
  delete [] A;
  delete [] Pos;
  for ( i = 0; i < Multiplicity; i++) delete [] Force[i];
  delete [] Force;
  delete [] ForceS;
  delete [] SavedVelo;

  return;
}


void G4StatMF::CalculateFragmentsMomentum(const G4int & INET,
					  const G4int & NFrags,
					  const G4double & T,
					  const G4double & TotKineticE,
					  G4ThreeVector * Momentum)
  // Calculates fragments momentum components at the breakup instant.
  // Fragment kinetic energies will be calculated according to the
  // Boltzamann distribution at given temperature.
{
  if (INET <= 0) return;
  G4int NFrags1 = NFrags;
  G4int NFragsM = NFrags + INET;

  if (INET == 1) {
    // INET == 1 only one fragment and absolute moment value will be
    Momentum[NFrags1] =
		   
IsotropicVector(sqrt(2.*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*
					     theSim->GetFragmentA(NFrags1)*TotKineticE));
  } else if (INET == 2) {
    Momentum[NFrags1] = IsotropicVector(sqrt(2.*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*
					     (theSim->GetFragmentA(NFrags1)*theSim->GetFragmentA(NFragsM))/
					     (theSim->GetFragmentA(NFrags1)+theSim->GetFragmentA(NFragsM))*
					     TotKineticE));
    Momentum[NFragsM] = Momentum[NFrags1];
  } else { // INET > 2
    // Sample kinetic energy sccording to Boltzmann distribution
    // Calculate the absolute values of fragments momenta
    // Calculate fragments momenta components (using isotropicl angular distrib.)
    // Sum fragments kinetic energies and fragments momentum components to check constraints.
    G4double EE = 0.0, Em = 0.0;
    G4int i1 = 0, i2 = 0;
    G4ThreeVector p;
    G4double Esum = 0.0;
    G4ThreeVector Psum(0.0,0.0,0.0);
    do {
      G4int NFragsM2 = NFragsM - 2;
      G4double FEMT = sqrt(0.5*T)*exp(-0.5);

      for (G4int i = NFrags1; i < NFragsM2; i++) {
	for (;;) {
	  G4double E = G4UniformRand()*9.0*T;
	  G4double FE = sqrt(E)*exp(-E/T);
	  G4double FErand = G4UniformRand()*FEMT;
	  if (FErand <= FE) {
	    Momentum[i] = IsotropicVector(sqrt(2.0*E*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*theSim->GetFragmentA(i)));
	    Esum += E;
	    Psum += Momentum[i];
	    break;
	  }
	}
      }
      // calculate momenta of two last fragments
      // to satisfy constraint
      i1 = NFragsM2;
      i2 = NFragsM2 + 1;
      p = -Psum;
      EE = TotKineticE - Esum;
      // Kinetic energy EE should be shared between two last fragments
      Em = p.mag2()/((G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*(theSim->GetFragmentA(i1)+theSim->GetFragmentA(i2)));
    } while (EE <= Em);
    G4double H = 1.0 + theSim->GetFragmentA(i2)/theSim->GetFragmentA(i1);
    G4double CTM12 = H*(1.0 - 2.0*theSim->GetFragmentA(i2)*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*EE/p.mag2());
    G4double CosTheta1;
    G4int zn;
    for (;;) {
      do {
	CosTheta1 = 1.0 - 2.0*G4UniformRand();
      } while (CosTheta1*CosTheta1 < CTM12);
      if (CTM12 < 0.0) {
	zn = 1.0;
	break;
      } else {
	if (CosTheta1 < 0.0) continue;
	else {
	  if (G4UniformRand() <= 0.5) {
	    zn = -1.0;
	    break;
	  } else {
	    zn = 1.0;
	    break;
	  }
	}
      }
    }
    G4double P1 = (p.mag()*CosTheta1+zn*sqrt(p.mag2()*CosTheta1*CosTheta1-p.mag2()*CTM12))/H;
    G4double P2 = sqrt(P1*P1+p.mag2() - 2.0*P1*p.mag()*CosTheta1);
    G4double Phi = twopi*G4UniformRand();
    G4double SinTheta1 = sqrt(1.0 - CosTheta1*CosTheta1);
    G4double CosPhi1 = cos(Phi);
    G4double SinPhi1 = sin(Phi);
    G4double CosPhi2 = -CosPhi1;
    G4double SinPhi2 = -SinPhi1;
    G4double CosTheta2 = (p.mag2() + P2*P2 - P1*P1)/(2.0*p.mag()*P2);
    G4double SinTheta2 = 0.0;
    if (CosTheta2 > -1.0 && CosTheta2 < 1.0) SinTheta2 = sqrt(1.0 - CosTheta2*CosTheta2);

    G4ThreeVector Pi(P1*SinTheta1*CosPhi1,P1*SinTheta1*SinPhi1,P1*CosTheta1);
    G4ThreeVector Pj(P2*SinTheta2*CosPhi2,P2*SinTheta2*SinPhi2,P2*CosTheta2);

    G4ThreeVector a = p;
    G4ThreeVector b(1.0,0.0,0.0);
    
    Momentum[i1] = Rotor(Pi,a,b);

    Momentum[i2] = Rotor(Pj,a,b);

    Psum += Momentum[i1] + Momentum[i2];
    Esum += Momentum[i1].mag2()/(2.0*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*theSim->GetFragmentA(i1))
      + Momentum[i2].mag2()/(2.0*(G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1)/MeV)*theSim->GetFragmentA(i2));
  }

  return;
}



G4ThreeVector G4StatMF::Rotor(const G4ThreeVector & P,
			      const G4ThreeVector & A,
			      const G4ThreeVector & B)
  // Rotates a 3-vector P to close momentum triangle P + A + B = 0
{
  G4double ScalarProd = A * B;
  G4double Alpha1 = ScalarProd/A.mag();
  G4double Alpha2 = sqrt(B.mag2() - Alpha1*Alpha1);
  G4ThreeVector NewV(A.y()*B.z()-A.z()*B.y(),
		     A.z()*B.x()-A.x()*B.z(),
		     A.x()*B.y()-A.y()*B.x());

  G4ThreeVector TheRotatedVector;

  TheRotatedVector.setX(P.x()*B.x()/Alpha2 +
			(P.z()-Alpha1*P.x()/Alpha2)*A.x()/A.mag() +
			(P.y()*A.x())/(Alpha2*A.mag()));


  TheRotatedVector.setX(P.x()*B.y()/Alpha2 +
			(P.z()-Alpha1*P.x()/Alpha2)*A.y()/A.mag() +
			(P.y()*A.y())/(Alpha2*A.mag()));

  TheRotatedVector.setX(P.x()*B.z()/Alpha2 +
			(P.z()-Alpha1*P.x()/Alpha2)*A.z()/A.mag() +
			(P.y()*A.z())/(Alpha2*A.mag()));

  return TheRotatedVector;
}




G4double G4StatMF::CalculateFragmentExcitationEnergy(const G4int & index, const G4double & T)
{

  G4double EvapEnergy = theSim->GetFragmentA(index)*T*T/theSim->GetFragmentInvLevelDensity(index);
  
  // For alpha particles
  if (theSim->GetFragmentA(index) == 4) return EvapEnergy;
  else {
    // Term connected with surface energy
    G4double ESurf;
    if (theSim->DBetaDT(T) == 0.0) ESurf = 0.0;
    else ESurf = (theSim->GetFragmentA(index)*T*T*2.5*G4StatMFParameters::GetBeta0())/
	   (G4StatMFParameters::GetCriticalTemp()*G4StatMFParameters::GetCriticalTemp()*
	    pow(theSim->GetFragmentA(index),1.0/3.0));
    EvapEnergy += ESurf;
  }
  return EvapEnergy;
}




G4ThreeVector G4StatMF::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*cos(Phi)*SinTheta,
		       Magnitude*cos(Phi)*CosTheta,
		       Magnitude*sin(Phi));
  return Vector;
}
