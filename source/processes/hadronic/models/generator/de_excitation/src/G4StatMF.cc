
#include "G4StatMF.hh"



// Default constructor
G4StatMF::G4StatMF() : _theEnsemble(0) {}


// Destructor
G4StatMF::~G4StatMF() {}; //{if (_theEnsemble != 0) delete _theEnsemble;}


// Copy constructor
G4StatMF::G4StatMF(const G4StatMF & right)
{
	G4Exception("G4StatMF::copy_constructor meant to not be accessable");
}


// Operators

G4StatMF & G4StatMF::operator=(const G4StatMF & right)
{
	G4Exception("G4StatMF::operator= meant to not be accessable");
	return *this;
}


G4bool G4StatMF::operator==(const G4StatMF & right)
{
	G4Exception("G4StatMF::operator== meant to not be accessable");
	return false;
}
 

G4bool G4StatMF::operator!=(const G4StatMF & right)
{
	G4Exception("G4StatMF::operator!= meant to not be accessable");
	return true;
}





G4FragmentVector * G4StatMF::BreakItUp(const G4Fragment &theFragment)
{
// 	G4FragmentVector * theResult = new G4FragmentVector;

	if (theFragment.GetExcitationEnergy() <= 0.0) {
		G4FragmentVector * theResult = new G4FragmentVector;
		theResult->insert(new G4Fragment(theFragment));
		return 0;
	}


	// Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
	// and M_0 = 3.3 for A <= 110
	G4double MaxAverageMultiplicity = 2.6;
	if (theFragment.GetA() <= 110) MaxAverageMultiplicity = 3.3;
	
	// We'll use two kinds of ensembles
	G4StatMFMicroCanonical * theMicrocanonicalEnsemble = 0;
  	G4StatMFMacroCanonical * theMacrocanonicalEnsemble = 0;

	
	//-------------------------------------------------------
	// Direct simulation part (Microcanonical ensemble)
	//-------------------------------------------------------
  
	// Microcanonical ensemble initialization 
	theMicrocanonicalEnsemble = new G4StatMFMicroCanonical(theFragment);

	G4int Iterations = 0;
	G4double Temperature = 0.0;
	G4double EnergyCoulomb = 0.0;
  
	G4bool FirstTime = true;
	G4StatMFChannel * theChannel = 0;

	G4bool ChannelOk;
	do {  // Try to de-excite as much as 10 times
		do {
		
			G4double theMeanMult = theMicrocanonicalEnsemble->GetMeanMultiplicity();
			if (theMeanMult <= MaxAverageMultiplicity) {
//				G4cout << "MICROCANONICAL" << G4endl;
				// choose fragments atomic numbers and charges from direct simulation
				theChannel = theMicrocanonicalEnsemble->ChooseAandZ(theFragment);
				_theEnsemble = theMicrocanonicalEnsemble;
			} else {
				//-------------------------------------------------
				// Non direct simulation part (Macrocanonical Ensemble)
				//-------------------------------------------------
				if (FirstTime) {
					// Macrocanonical ensemble initialization 
					theMacrocanonicalEnsemble = new G4StatMFMacroCanonical(theFragment);
					_theEnsemble = theMacrocanonicalEnsemble;
					FirstTime = false;
				}
//				G4cout << "MACROCANONICAL" << G4endl;
				// Select calculated fragment total multiplicity, 
				// fragment atomic numbers and fragment charges.
				theChannel = theMacrocanonicalEnsemble->ChooseAandZ(theFragment);
			}
      	
			if (!(ChannelOk = theChannel->CheckFragments())) delete theChannel; 
			
		} while (!ChannelOk);


		if (theChannel->GetMultiplicity() <= 1) {
			G4FragmentVector * theResult = new G4FragmentVector;
			theResult->insert(new G4Fragment(theFragment));
			delete theMicrocanonicalEnsemble;
			if (theMacrocanonicalEnsemble != 0) delete theMacrocanonicalEnsemble;
			delete theChannel;
			return theResult;
		}

		//--------------------------------------
		// Second part of simulation procedure.
		//--------------------------------------

		// Find temperature of breaking channel.
		Temperature = _theEnsemble->GetMeanTemperature(); // Initial value for Temperature 

		if (FindTemperatureOfBreakingChannel(theFragment,theChannel,Temperature)) break;

	} while (Iterations++ < 10);


	// If Iterations >= 10 means that we couldn't solve for temperature
	if (Iterations >= 10) G4Exception("G4StatMF::BreakItUp: Was not possible to solve for temperature of breaking channel");


	G4FragmentVector * theResult = theChannel->GetFragments(theFragment.GetA(),theFragment.GetZ(),Temperature);

	
	
	// ~~~~~~ Energy conservation Patch !!!!!!!!!!!!!!!!!!!!!!
	// Original nucleus 4-momentum in CM system
	G4LorentzVector InitialMomentum(theFragment.GetMomentum());
	InitialMomentum.boost(-InitialMomentum.boostVector());
	G4double ScaleFactor;
	do {
		G4double FragmentsEnergy = 0.0;
		G4int j = 0;
		for (j = 0; j < theResult->entries(); j++) FragmentsEnergy += theResult->at(j)->GetMomentum().e();
		ScaleFactor = InitialMomentum.e()/FragmentsEnergy;
		G4ThreeVector ScaledMomentum(0.0,0.0,0.0);
		for (j = 0; j < theResult->entries(); j++) {
			ScaledMomentum = ScaleFactor*theResult->at(j)->GetMomentum().vect();
			G4double Mass = theResult->at(j)->GetMomentum().m();
			G4LorentzVector NewMomentum;
			NewMomentum.setVect(ScaledMomentum);
			NewMomentum.setE(sqrt(ScaledMomentum.mag2()+Mass*Mass));
			theResult->at(j)->SetMomentum(NewMomentum);		
		}
	} while (ScaleFactor > 1.0+1.e-10);
	// ~~~~~~ End of patch !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	// Perform Lorentz boost
	for (G4int i = 0; i < theResult->entries(); i++) {
		G4LorentzVector FourMom = theResult->at(i)->GetMomentum();
		FourMom.boost(theFragment.GetMomentum().boostVector());
		theResult->at(i)->SetMomentum(FourMom);
	}

	// garbage collection
	delete theMicrocanonicalEnsemble;
	if (theMacrocanonicalEnsemble != 0) delete theMacrocanonicalEnsemble;
	delete theChannel;

	return theResult;
}


G4bool G4StatMF::FindTemperatureOfBreakingChannel(const G4Fragment & theFragment,
						  const G4StatMFChannel * aChannel,
						  G4double & Temperature)
  // This finds temperature of breaking channel.
{
	G4double A = theFragment.GetA();
	G4double Z = theFragment.GetZ();
	G4double U = theFragment.GetExcitationEnergy();

	G4double T = G4std::max(Temperature,0.0012*MeV);

	G4double Ta = T;
	G4double Tb = T;
        
   
	G4double TotalEnergy = CalcEnergy(A,Z,aChannel,T);

	G4double Da = (U - TotalEnergy)/U;
	G4double Db = 0.0;
    
	// bracketing the solution
	if (Da == 0.0) {
		Temperature = T;
		return true;
	} else if (Da < 0.0) {
		do {
			Tb -= 0.5 * abs(Tb);
			T = Tb;
			if (Tb < 0.001*MeV) return false;
            
			TotalEnergy = CalcEnergy(A,Z,aChannel,T);
        
			Db = (U - TotalEnergy)/U;
		} while (Db < 0.0);

	} else {
		do {
			Tb += 0.5 * abs(Tb);
			T = Tb;
        
			TotalEnergy = CalcEnergy(A,Z,aChannel,T);
        
			Db = (U - TotalEnergy)/U;
		} while (Db > 0.0); 
	}
    
	G4double eps = 1.0e-14 * abs(Tb-Ta);
	//G4double eps = 1.0e-3 ;
    
	// Start the bisection method
	for (G4int j = 0; j < 1000; j++) {
		G4double Tc =  (Ta+Tb)/2.0;
		if (abs(Ta-Tc) <= eps) {
			Temperature = Tc;
			return true;
		}
        
		T = Tc;
        
		TotalEnergy = CalcEnergy(A,Z,aChannel,T);
        
		G4double Dc = (U - TotalEnergy)/U; 
        
		if (Dc == 0.0) {
			Temperature  = Tc;
			return true;
		}
        
		if (Da*Dc < 0.0) {
			Tb = Tc;
			Db = Dc;
		} else {
			Ta = Tc;
			Da = Dc;
		}
	}
    
	Temperature  = (Ta+Tb)/2.0;
	return false;
// ***********************************************************
//   const G4double HT = 0.5;
//   G4double H = 0.0;
//   G4int counter = 0;
//   G4int id = 0;
//   G4int K1 = 0, K2 = 0;

//   do {

//     G4double TotalEnergy = CalcEnergy(A,Z,aChannel,T);
// 
//     G4double D = (U - TotalEnergy)/U;
//     if (abs(D) < 0.003) {
//       Temperature = T;
//       return true;
//     }
//     counter++;
// 
//     if (D < 0.0) H = -HT;
//     else H = HT;
// 
//     if (D <= 0.0) {
//       for (;;) {
// 	    K1 = 1;
// 	    if (K1 ==1 && K2 == 1) id++;
// 	    if (id > 30) return false;
// 	    T += H/pow(2.0,id);
// 	    if ( T >= 0.001) break;
// 	    K2 = 1;
// 	    H = HT;
//       }
//     } else {
//       for (;;) {
// 	    K2 = 1;
// 	    if (K1 ==1 && K2 == 1) id++;
// 	    if (id > 30) return false;
// 	    T += H/pow(2.0,id);
// 	    if ( T >= 0.001) break;
// 	    K1 = 1;
// 	    H = HT;
//       }
//     }
//   } while (counter <= 120);
//   return false;
}



G4double G4StatMF::CalcEnergy(const G4double A, const G4double Z, const G4StatMFChannel * aChannel,
										const G4double T)
{
	G4double MassExcess0 = G4NucleiProperties::GetMassExcess(A,Z);
	
	G4double Coulomb = (3./5.)*(elm_coupling*Z*Z)*pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.)/
	 						(G4StatMFParameters::Getr0()*pow(A,1./3.));

   G4double ChannelEnergy = aChannel->GetFragmentsEnergy(T);
	
	return -MassExcess0 + Coulomb + ChannelEnergy;

}



