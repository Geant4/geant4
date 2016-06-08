#include "G4StatMFChannel.hh"



// Copy constructor
G4StatMFChannel::G4StatMFChannel(const G4StatMFChannel & right)
{
	G4Exception("G4StatMFChannel::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFChannel & G4StatMFChannel::
operator=(const G4StatMFChannel & right)
{
	G4Exception("G4StatMFChannel::operator= meant to not be accessable");
	return *this;
}


G4bool G4StatMFChannel::operator==(const G4StatMFChannel & right) const
{
//	G4Exception("G4StatMFChannel::operator== meant to not be accessable");
	return false;
}
 

G4bool G4StatMFChannel::operator!=(const G4StatMFChannel & right) const 
{
//	G4Exception("G4StatMFChannel::operator!= meant to not be accessable");
	return true;
}


G4bool G4StatMFChannel::CheckFragments(void)
{
	for (G4int i = 0; i < _theFragments.entries(); i++) {
		G4int A = G4int(_theFragments(i)->GetA());
		G4int Z = G4int(_theFragments(i)->GetZ());
		if (A > 1 && (Z >= A || Z <= 0) || (A==1 && Z > A) || A <= 0) return false;
	}
	
	return true;
}




void G4StatMFChannel::CreateFragment(const G4double A, const G4double Z)
// Create a new fragment.
// Fragments are automatically sorted: first charged fragments, 
// then neutral ones.
{
	if (Z <= 0.5) {
		_theFragments.insert(new G4StatMFFragment(A,Z));
		_NumOfNeutralFragments++;
	} else {
		_theFragments.prepend(new G4StatMFFragment(A,Z));
		_NumOfChargedFragments++;
	}
	
	return;
}


G4double G4StatMFChannel::GetFragmentsCoulombEnergy(void)
{
	G4double Coulomb = 0.0;
	for (G4int i = 0;i < _theFragments.entries(); i++)
						Coulomb += _theFragments(i)->GetCoulombEnergy();
	return Coulomb;
}



G4double G4StatMFChannel::GetFragmentsEnergy(const G4double T) const
{
	G4double Energy = 0.0;
	
	G4double TranslationalEnergy = (3./2.)*T*_theFragments.entries();
	
	for (G4int i = 0;i < _theFragments.entries(); i++)
						Energy += _theFragments(i)->GetEnergy(T);
	return Energy + TranslationalEnergy;	
}

G4FragmentVector * G4StatMFChannel::GetFragments(const G4double anA, 
																const G4double anZ,
																const G4double T)
	//
{

	// calculate momenta of charged fragments  
	CoulombImpulse(anA,anZ,T);
	
	// calculate momenta of neutral fragments
	FragmentsMomenta(_NumOfNeutralFragments, _NumOfChargedFragments, T);


	G4FragmentVector * theResult = new G4FragmentVector;
	for (G4int i = 0; i < _theFragments.entries(); i++)
			theResult->insert(_theFragments(i)->GetFragment(T));

	return theResult;

}



void G4StatMFChannel::CoulombImpulse(const G4double anA, const G4double anZ, const G4double T)
	// Aafter breakup, fragments fly away under Coulomb field.
	// This method calculates asymptotic fragments momenta.
{
	// First, we have to place the fragments inside of the original nucleus volume
	PlaceFragments(anA);
	
	// Second, we sample initial charged fragments momenta. There are
	// _NumOfChargedFragments charged fragments and they start at the begining
	// of the vector _theFragments (i.e. 0) 
	FragmentsMomenta(_NumOfChargedFragments, 0, T);

	// Third, we have to figure out the asymptotic momenta of charged fragments 
	// For taht we have to solve equations of motion for fragments
	SolveEqOfMotion(anA,anZ,T);

	return;
}



void G4StatMFChannel::PlaceFragments(const G4double anA)
	// This gives the position of fragments at the breakup instant. 
	// Fragments positions are sampled inside prolongated ellipsoid.
{
	const G4double R0 = G4StatMFParameters::Getr0();
	const G4double Rsys = 2.0*R0*pow(anA,1./3.);

	G4bool TooMuchIterations;
	do {
		TooMuchIterations = false;

		// Sample the position of the first fragment
		G4double R = (Rsys - R0*pow(_theFragments(0)->GetA(),1./3.))*
							pow(G4UniformRand(),1./3.);
		_theFragments(0)->SetPosition(IsotropicVector(R));
		

		// Sample the position of the remaining fragments
		G4bool ThereAreOverlaps = false;
		for (G4int i = 1; i < _theFragments.entries(); i++) {
			G4int counter = 0;
			do {
				R = (Rsys - R0*pow(_theFragments(i)->GetA(),1./3.))*
						pow(G4UniformRand(),1./3.);
				_theFragments(i)->SetPosition(IsotropicVector(R));
			
				// Check that there are not overlapping fragments
				for (G4int j = 0; j < i; j++) {
					G4ThreeVector FragToFragVector = 
						_theFragments(i)->GetPosition() - _theFragments(j)->GetPosition();
					
					G4double Rmin = R0*(pow(_theFragments(i)->GetA(),1./3.) +
										  pow(_theFragments(j)->GetA(),1./3));
					if (ThereAreOverlaps = (FragToFragVector.mag2() < Rmin*Rmin)) break;
				}
				counter++;
			} while (ThereAreOverlaps && counter < 1000);
			
			if (counter >= 1000) {
				TooMuchIterations = true;
				break;
			} 
		}
	} while (TooMuchIterations);
	
	return;
}


void G4StatMFChannel::FragmentsMomenta(const G4int NF, const G4int idx,
													const G4double T)
// Calculate fragments momenta at the breakup instant
// Fragment kinetic energies are calculated according to the 
// Boltzmann distribution at given temperature.
// NF is number of fragments
// idx is index of first fragment
{
	G4double KinE = (3./2.)*T*G4double(NF);
	
	G4ThreeVector p;
	
	if (NF <= 0) return;
	else if (NF == 1) {
		// We have only one fragment to deal with
		p = IsotropicVector(sqrt(2.0*_theFragments(idx)->GetNuclearMass()*KinE));
		_theFragments(idx)->SetMomentum(p);
	} else if (NF == 2) {
		// We have only two fragment to deal with
		G4double M1 = _theFragments(idx)->GetNuclearMass();
		G4double M2 = _theFragments(idx+1)->GetNuclearMass();
		p = IsotropicVector(sqrt(2.0*KinE*(M1*M2)/(M1+M2)));		
		_theFragments(idx)->SetMomentum(p);
		_theFragments(idx+1)->SetMomentum(-p);
	} else {
		// We have more than two fragments
		G4double AvailableE;
		G4int i1,i2;
		G4double SummedE;
		G4ThreeVector SummedP;
		do {
			// Fisrt sample momenta of NF-2 fragments 
			// according to Boltzmann distribution
			AvailableE = 0.0;
			SummedE = 0.0;
			SummedP.setX(0.0);SummedP.setY(0.0);SummedP.setZ(0.0);
			for (G4int i = idx; i < idx+NF-2; i++) {
				G4double E;
				G4double RandE;
				G4double Boltzmann;
				do {
					E = 9.0*T*G4UniformRand();
					Boltzmann = sqrt(E)*exp(-E/T);
					RandE = sqrt(T/2.)*exp(-0.5)*G4UniformRand();
				} while (RandE > Boltzmann);
				p = IsotropicVector(sqrt(2.0*E*_theFragments(i)->GetNuclearMass()));
				_theFragments(i)->SetMomentum(p);
				SummedE += E;
				SummedP += p;
			}	
			// Calculate momenta of last two fragments in such a way
			// that constraints are satisfied
			i1 = idx+NF-2;  // before last fragment index
			i2 = idx+NF-1;  // last fragment index
			p = -SummedP;
			AvailableE = KinE - SummedE;
			// Available Kinetic Energy should be shared between two last fragments
		} while (AvailableE <= p.mag2()/(2.0*(_theFragments(i1)->GetNuclearMass()+
														  _theFragments(i2)->GetNuclearMass())));		
		
//		G4double H = 1.0 + _theFragments(i2)->GetA()/_theFragments(i1)->GetA();
		G4double H = 1.0 + _theFragments(i2)->GetNuclearMass()/_theFragments(i1)->GetNuclearMass();
		G4double CTM12 = H*(1.0 - 2.0*_theFragments(i2)->GetNuclearMass()*AvailableE/p.mag2());
		G4double CosTheta1;
		G4double Sign;
		do {
			do {
				CosTheta1 = 1.0 - 2.0*G4UniformRand();
			} while (CosTheta1*CosTheta1 < CTM12);
		} while (CTM12 >= 0.0 && CosTheta1 < 0.0);

		if (CTM12 < 0.0) Sign = 1.0;
		else if (G4UniformRand() <= 0.5) Sign = -1.0;
		else Sign = 1.0;
		
		
		G4double P1 = (p.mag()*CosTheta1+Sign*sqrt(p.mag2()*(CosTheta1*CosTheta1-CTM12)))/H;
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
	
		G4ThreeVector p1(P1*SinTheta1*CosPhi1,P1*SinTheta1*SinPhi1,P1*CosTheta1);
		G4ThreeVector p2(P2*SinTheta2*CosPhi2,P2*SinTheta2*SinPhi2,P2*CosTheta2);
		G4ThreeVector b(1.0,0.0,0.0);
		
		p1 = RotateMomentum(p,b,p1);
		p2 = RotateMomentum(p,b,p2);
	
		SummedP += p1 + p2;
		SummedE += p1.mag2()/(2.0*_theFragments(i1)->GetNuclearMass()) + 
					  p2.mag2()/(2.0*_theFragments(i2)->GetNuclearMass());		
		
		_theFragments(i1)->SetMomentum(p1);
		_theFragments(i2)->SetMomentum(p2);
		
	}

	return;
}


void G4StatMFChannel::SolveEqOfMotion(const G4double anA, const G4double anZ, const G4double T)
// This method will find a solution of Newton's equation of motion
// for fragments in the self-consistent time-dependent Coulomb field
{
	G4double CoulombEnergy = (3./5.)*(elm_coupling*anZ*anZ)*
							pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.)/
	 						(G4StatMFParameters::Getr0()*pow(anA,1./3.))
							- GetFragmentsCoulombEnergy();
	if (CoulombEnergy <= 0.0) return;

	G4int Iterations = 0;
	G4double TimeN = 0.0;
	G4double TimeS = 0.0;
	G4double DeltaTime = 10.0;

	G4ThreeVector * Pos = new G4ThreeVector[_NumOfChargedFragments];
	G4ThreeVector * Vel = new G4ThreeVector[_NumOfChargedFragments];
	G4ThreeVector * Accel = new G4ThreeVector[_NumOfChargedFragments];

	G4int i;
	for (i = 0; i < _NumOfChargedFragments; i++) {
		Vel[i] = (1.0/(_theFragments(i)->GetNuclearMass()))*
						_theFragments(i)->GetMomentum();
		Pos[i] = _theFragments(i)->GetPosition();
	}
  
	do {

		G4ThreeVector distance;
		G4ThreeVector force;

		for (i = 0; i < _NumOfChargedFragments; i++) {
			force.setX(0.0); force.setY(0.0); force.setZ(0.0);
			for (G4int j = 0; j < _NumOfChargedFragments; j++) {
				if (i != j) {
					distance = Pos[i] - Pos[j];
					force += (elm_coupling*(_theFragments(i)->GetZ()*_theFragments(j)->GetZ())/
										(distance.mag2()*distance.mag()))*distance;
				}
			}
			Accel[i] = (1./(_theFragments(i)->GetNuclearMass()))*force;
		}

		TimeN = TimeS + DeltaTime;
	
		G4ThreeVector SavedVel;
		for ( i = 0; i < _NumOfChargedFragments; i++) {
			SavedVel = Vel[i];
			Vel[i] += Accel[i]*(TimeN-TimeS);
			Pos[i] += (SavedVel+Vel[i])*(TimeN-TimeS)*0.5;
		}

// 		if (Iterations >= 50 && Iterations < 75) DeltaTime = 4.;
// 		else if (Iterations >= 75) DeltaTime = 10.;

		TimeS = TimeN;

	} while (Iterations++ < 100);
	
	// Summed fragment kinetic energy
	G4double TotalKineticEnergy = 0.0;
	for (i = 0; i < _NumOfChargedFragments; i++) TotalKineticEnergy += 
		_theFragments(i)->GetNuclearMass()*0.5*Vel[i].mag2();
	
	// Scaling of fragment velocities
	G4double KineticEnergy = (3./2.)*G4double(_theFragments.entries())*T;
	G4double Eta = ( CoulombEnergy + KineticEnergy ) / TotalKineticEnergy;
	for (i = 0; i < _NumOfChargedFragments; i++) Vel[i] *= Eta;
	
	
	// Finally calculate fragments momenta
	for (i = 0; i < _NumOfChargedFragments; i++) 
				_theFragments(i)->SetMomentum(_theFragments(i)->GetNuclearMass()*Vel[i]);

	
	// garbage collection
	delete [] Pos;
	delete [] Vel;
	delete [] Accel;

	return;
}



G4ThreeVector G4StatMFChannel::RotateMomentum(G4ThreeVector Pa,
						      G4ThreeVector V, G4ThreeVector P)
// Rotates a 3-vector P to close momentum triangle Pa + V + P = 0
{
  G4ThreeVector U = Pa.unit();
  
  G4double Alpha1 = U * V;
  
  G4double Alpha2 = sqrt(V.mag2() - Alpha1*Alpha1);

  G4ThreeVector N = (1./Alpha2)*U.cross(V);

  G4ThreeVector RotatedMomentum(
				     ( (V.x() - Alpha1*U.x())/Alpha2 ) * P.x() + N.x() * P.y() + U.x() * P.z(),
				     ( (V.y() - Alpha1*U.y())/Alpha2 ) * P.x() + N.y() * P.y() + U.y() * P.z(),
				     ( (V.z() - Alpha1*U.z())/Alpha2 ) * P.x() + N.z() * P.y() + U.z() * P.z()
				     );
  return RotatedMomentum;
}





G4ThreeVector G4StatMFChannel::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vector with a magnitud given by Magnitude.
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
