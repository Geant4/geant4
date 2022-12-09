//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMF.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4PhysicsModelCatalog.hh"

G4StatMF::G4StatMF()
{
  _secID = G4PhysicsModelCatalog::GetModelID("model_G4StatMF");
}

G4StatMF::~G4StatMF() {}

G4FragmentVector* G4StatMF::BreakItUp(const G4Fragment &theFragment)
{
  if (theFragment.GetExcitationEnergy() <= 0.0) {
    return nullptr;
  }

  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
  // and M_0 = 3.3 for A <= 110
  G4double MaxAverageMultiplicity = 
    G4StatMFParameters::GetMaxAverageMultiplicity(theFragment.GetA_asInt());

	
    // We'll use two kinds of ensembles
  G4StatMFMicroCanonical * theMicrocanonicalEnsemble = 0;
  G4StatMFMacroCanonical * theMacrocanonicalEnsemble = 0;
	
  //-------------------------------------------------------
  // Direct simulation part (Microcanonical ensemble)
  //-------------------------------------------------------
  
  // Microcanonical ensemble initialization 
  theMicrocanonicalEnsemble = new G4StatMFMicroCanonical(theFragment);

  G4int Iterations = 0;
  G4int IterationsLimit = 100000;
  G4double Temperature = 0.0;
  
  G4bool FirstTime = true;
  G4StatMFChannel * theChannel = 0;
 
  G4bool ChannelOk;
  do {  // Try to de-excite as much as IterationLimit permits
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
      
      ChannelOk = theChannel->CheckFragments();
      if (!ChannelOk) delete theChannel; 
      
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
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
    Temperature = _theEnsemble->GetMeanTemperature(); // Initial guess for Temperature 
 
    if (FindTemperatureOfBreakingChannel(theFragment,theChannel,Temperature)) break;
 
    // Do not forget to delete this unusable channel, for which we failed to find the temperature,
    // otherwise for very proton-reach nuclei it would lead to memory leak due to large 
    // number of iterations. N.B. "theChannel" is created in G4StatMFMacroCanonical::ChooseZ()

    // G4cout << " Iteration # " << Iterations << " Mean Temperature = " << Temperature << G4endl;    

    delete theChannel;    

    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while (Iterations++ < IterationsLimit );

  // If Iterations >= IterationsLimit means that we couldn't solve for temperature
  if (Iterations >= IterationsLimit) 
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMF::BreakItUp: Was not possible to solve for temperature of breaking channel");

  G4FragmentVector * theResult = theChannel->
    GetFragments(theFragment.GetA_asInt(),theFragment.GetZ_asInt(),Temperature);
  	
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
      NewMomentum.setE(std::sqrt(ScaledMomentum.mag2()+Mass*Mass));
      (*j)->SetMomentum(NewMomentum);		
    }
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  } while (ScaleFactor > 1.0+1.e-5 && std::abs(ScaleFactor-SavedScaleFactor)/ScaleFactor > 1.e-10);
  // ~~~~~~ End of patch !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  // Perform Lorentz boost
  G4FragmentVector::iterator i;
  for (i = theResult->begin(); i != theResult->end(); i++) {
    G4LorentzVector FourMom = (*i)->GetMomentum();
    FourMom.boost(theFragment.GetMomentum().boostVector());
    (*i)->SetMomentum(FourMom);
    (*i)->SetCreatorModelID(_secID);
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
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  G4double U = theFragment.GetExcitationEnergy();
  
  G4double T = std::max(Temperature,0.0012*MeV);  
  G4double Ta = T;
  G4double TotalEnergy = CalcEnergy(A,Z,aChannel,T);
  
  G4double Da = (U - TotalEnergy)/U;
  G4double Db = 0.0;
  
  // bracketing the solution
  if (Da == 0.0) {
    Temperature = T;
    return true;
  } else if (Da < 0.0) {
    do {
      T *= 0.5;
      if (T < 0.001*MeV) return false;
      
      TotalEnergy = CalcEnergy(A,Z,aChannel,T);
      
      Db = (U - TotalEnergy)/U;
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (Db < 0.0);
    
  } else {
    do {
      T *= 1.5;
      
      TotalEnergy = CalcEnergy(A,Z,aChannel,T);
      
      Db = (U - TotalEnergy)/U;
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (Db > 0.0); 
  }
  
  G4double eps = 1.0e-14 * std::abs(T-Ta);
  //G4double eps = 1.0e-3 ;
  
  // Start the bisection method
  for (G4int j = 0; j < 1000; j++) {
    G4double Tc =  (Ta+T)*0.5;
    if (std::abs(Ta-Tc) <= eps) {
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
      T  = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    }
  }
  
  Temperature  = (Ta+T)*0.5;
  return false;
}

G4double G4StatMF::CalcEnergy(G4int A, G4int Z, const G4StatMFChannel * aChannel,
			      G4double T)
{
  G4double MassExcess0 = G4NucleiProperties::GetMassExcess(A,Z);
  G4double ChannelEnergy = aChannel->GetFragmentsEnergy(T);
  return -MassExcess0 + G4StatMFParameters::GetCoulomb() + ChannelEnergy;
}



