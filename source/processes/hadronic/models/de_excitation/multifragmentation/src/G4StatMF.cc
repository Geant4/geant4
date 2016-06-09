//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMF.cc,v 1.2 2003/11/03 17:53:05 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMF.hh"



// Default constructor
G4StatMF::G4StatMF() : _theEnsemble(0) {}


// Destructor
G4StatMF::~G4StatMF() {} //{if (_theEnsemble != 0) delete _theEnsemble;}


// Copy constructor
G4StatMF::G4StatMF(const G4StatMF & ) : G4VMultiFragmentation()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::copy_constructor meant to not be accessable");
}


// Operators

G4StatMF & G4StatMF::operator=(const G4StatMF & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMF::operator==(const G4StatMF & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMF::operator!=(const G4StatMF & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::operator!= meant to not be accessable");
    return true;
}





G4FragmentVector * G4StatMF::BreakItUp(const G4Fragment &theFragment)
{
  // 	G4FragmentVector * theResult = new G4FragmentVector;

  if (theFragment.GetExcitationEnergy() <= 0.0) {
    G4FragmentVector * theResult = new G4FragmentVector;
    theResult->push_back(new G4Fragment(theFragment));
    return 0;
  }


  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
  // and M_0 = 3.3 for A <= 110
  G4double MaxAverageMultiplicity = 
    G4StatMFParameters::GetMaxAverageMultiplicity(static_cast<G4int>(theFragment.GetA()));

	
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
  
  G4bool FirstTime = true;
  G4StatMFChannel * theChannel = 0;

  G4bool ChannelOk;
  do {  // Try to de-excite as much as 10 times
    do {
      
      G4double theMeanMult = theMicrocanonicalEnsemble->GetMeanMultiplicity();
      if (theMeanMult <= MaxAverageMultiplicity) {
	// G4cout << "MICROCANONICAL" << G4endl;
	// Choose fragments atomic numbers and charges from direct simulation
	theChannel = theMicrocanonicalEnsemble->ChooseAandZ(theFragment);
	_theEnsemble = theMicrocanonicalEnsemble;
      } else {
	//-----------------------------------------------------
	// Non direct simulation part (Macrocanonical Ensemble)
	//-----------------------------------------------------
	if (FirstTime) {
	  // Macrocanonical ensemble initialization 
	  theMacrocanonicalEnsemble = new G4StatMFMacroCanonical(theFragment);
	  _theEnsemble = theMacrocanonicalEnsemble;
	  FirstTime = false;
	}
	// G4cout << "MACROCANONICAL" << G4endl;
	// Select calculated fragment total multiplicity, 
	// fragment atomic numbers and fragment charges.
	theChannel = theMacrocanonicalEnsemble->ChooseAandZ(theFragment);
      }
      
      if (!(ChannelOk = theChannel->CheckFragments())) delete theChannel; 
      
    } while (!ChannelOk);
    
    
    if (theChannel->GetMultiplicity() <= 1) {
      G4FragmentVector * theResult = new G4FragmentVector;
      theResult->push_back(new G4Fragment(theFragment));
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
  if (Iterations >= 10) 
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::BreakItUp: Was not possible to solve for temperature of breaking channel");
  
  
  G4FragmentVector * theResult = theChannel->
    GetFragments(theFragment.GetA(),theFragment.GetZ(),Temperature);
  
  
	
  // ~~~~~~ Energy conservation Patch !!!!!!!!!!!!!!!!!!!!!!
  // Original nucleus 4-momentum in CM system
  G4LorentzVector InitialMomentum(theFragment.GetMomentum());
  InitialMomentum.boost(-InitialMomentum.boostVector());
  G4double ScaleFactor = 0.0;
  G4double SavedScaleFactor = 0.0;
  do {
    G4double FragmentsEnergy = 0.0;
    G4FragmentVector::iterator j;
    for (j = theResult->begin(); j != theResult->end(); j++) 
      FragmentsEnergy += (*j)->GetMomentum().e();
    SavedScaleFactor = ScaleFactor;
    ScaleFactor = InitialMomentum.e()/FragmentsEnergy;
    G4ThreeVector ScaledMomentum(0.0,0.0,0.0);
    for (j = theResult->begin(); j != theResult->end(); j++) {
      ScaledMomentum = ScaleFactor * (*j)->GetMomentum().vect();
      G4double Mass = (*j)->GetMomentum().m();
      G4LorentzVector NewMomentum;
      NewMomentum.setVect(ScaledMomentum);
      NewMomentum.setE(sqrt(ScaledMomentum.mag2()+Mass*Mass));
      (*j)->SetMomentum(NewMomentum);		
    }
  } while (ScaleFactor > 1.0+1.e-5 && abs(ScaleFactor-SavedScaleFactor)/ScaleFactor > 1.e-10);
  // ~~~~~~ End of patch !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  // Perform Lorentz boost
  G4FragmentVector::iterator i;
  for (i = theResult->begin(); i != theResult->end(); i++) {
    G4LorentzVector FourMom = (*i)->GetMomentum();
    FourMom.boost(theFragment.GetMomentum().boostVector());
    (*i)->SetMomentum(FourMom);
#ifdef PRECOMPOUND_TEST
    (*i)->SetCreatorModel(G4String("G4StatMF"));
#endif
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
  
  G4double T = std::max(Temperature,0.0012*MeV);
  
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
}



G4double G4StatMF::CalcEnergy(const G4double A, const G4double Z, const G4StatMFChannel * aChannel,
			      const G4double T)
{
    G4double MassExcess0 = G4NucleiProperties::GetMassExcess(static_cast<G4int>(A),static_cast<G4int>(Z));
	
    G4double Coulomb = (3./5.)*(elm_coupling*Z*Z)*pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.)/
	(G4StatMFParameters::Getr0()*pow(A,1./3.));

    G4double ChannelEnergy = aChannel->GetFragmentsEnergy(T);
	
    return -MassExcess0 + Coulomb + ChannelEnergy;

}



