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
// Multi-fragmentation 
// by V. Lara
//

#include "G4StatMF.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4PhysicsModelCatalog.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

G4StatMF::G4StatMF()
{
  theMicrocanonicalEnsemble = new G4StatMFMicroCanonical();
  theMacrocanonicalEnsemble = new G4StatMFMacroCanonical();
  //fSecID = G4PhysicsModelCatalog::GetModelID("model_G4StatMF");
}

G4StatMF::~G4StatMF()
{
  delete theMicrocanonicalEnsemble;
  delete theMacrocanonicalEnsemble;
}

G4FragmentVector* G4StatMF::BreakItUp(const G4Fragment& theFragment)
{
  // Maximun average multiplicity: M_0 = 2.6 for A ~ 200 
  // and M_0 = 3.3 for A <= 110
  G4int A = theFragment.GetA_asInt();
  G4double MaxAverageMultiplicity = G4StatMFParameters::GetMaxAverageMultiplicity(A);
  
  // Microcanonical ensemble - direct simulation
  theMicrocanonicalEnsemble->Initialise(theFragment);
  
  const G4int iLimit = 20;
  G4double Temperature = 0.0;
  
   G4StatMFChannel* theChannel = nullptr;
 
  for (G4int i=0; i<iLimit; ++i) {
    for (G4int j=0; j<iLimit; ++j) {
      
      G4double theMeanMult = theMicrocanonicalEnsemble->GetMeanMultiplicity();
      if (theMeanMult <= MaxAverageMultiplicity) {
	//G4cout << "MICROCANONICAL Nmean=" << theMeanMult
	//       << " i=" << i << " j=" << j << G4endl;
	// Choose fragments atomic numbers and charges from direct simulation
	theChannel = theMicrocanonicalEnsemble->ChooseAandZ(theFragment);
	fEnsemble = theMicrocanonicalEnsemble;
      } else {
	//-----------------------------------------------------
	// Non direct simulation part (Macrocanonical Ensemble)
	//-----------------------------------------------------
	// Macrocanonical ensemble initialization 
	theMacrocanonicalEnsemble->Initialise(theFragment);
	fEnsemble = theMacrocanonicalEnsemble;
	//G4cout << "MACROCANONICAL Nmean=" << theMeanMult
	//       << " i=" << i << " j=" << j << G4endl;
	// Select calculated fragment total multiplicity, 
	// fragment atomic numbers and fragment charges.
	theChannel = theMacrocanonicalEnsemble->ChooseAandZ(theFragment);
      }
      
      if (theChannel->CheckFragments()) { break; }
      delete theChannel; 
      theChannel = nullptr;
    }
        
    if (nullptr == theChannel || theChannel->GetMultiplicity() <= 1) {
      delete theChannel;
      theChannel = nullptr;
      break;
    }

    // G4cout << "   multiplicity=" << theChannel->GetMultiplicity() << G4endl;
    //--------------------------------------
    // Second part of simulation procedure.
    //--------------------------------------
    
    // Find temperature of breaking channel.
    Temperature = fEnsemble->GetMeanTemperature(); // Initial guess for Temperature 
 
    if (FindTemperatureOfBreakingChannel(theFragment, theChannel, Temperature)) {
      break;
    }
 
    // Do not forget to delete this unusable channel, for which we failed to find the temperature,
    // otherwise for very proton-reach nuclei it would lead to memory leak due to large 
    // number of iterations. N.B. "theChannel" is created in G4StatMFMacroCanonical::ChooseZ()

    // G4cout << " Iteration # " << Iterations << " Mean Temperature = " << Temperature << G4endl;
    delete theChannel;
    theChannel = nullptr;
  }

  // primary
  G4FragmentVector* theResult = nullptr;
  
  // no multi-fragmentation
  if (nullptr == theChannel || theChannel->GetMultiplicity() <= 1) {
    theResult = new G4FragmentVector();
    theResult->push_back(new G4Fragment(theFragment));
    delete theChannel;
    return theResult;
  }
  
  G4int Z = theFragment.GetZ_asInt();
  G4double m0 = theFragment.GetGroundStateMass() + theFragment.GetExcitationEnergy();
  auto bs = theFragment.GetMomentum().boostVector();
  	
  G4double etot = 0.0;
  G4double ekin = 0.0;
  for (G4int i=0; i<iLimit; ++i) {
    theResult = theChannel->GetFragments(A, Z, Temperature);
    if (nullptr == theResult) { continue; }
    etot = 0.0;
    ekin = 0.0;
    for (auto const & ptr : *theResult) {
      G4double e = ptr->GetMomentum().e();
      G4double m1 = ptr->GetGroundStateMass() + ptr->GetExcitationEnergy();
      etot += e;
      ekin += std::max(e - m1, 0.0);
    }
    // correction possible
    if (etot - m0 + ekin > 0.0 && ekin > 0.0) { break; }

    // new attemt required
    for (auto const & ptr : *theResult) {
      delete ptr;
    }
    delete theResult;
    theResult = nullptr;
  }
  delete theChannel;
  
  // no multi-fragmentation
  if (nullptr == theResult || ekin <= 0.0) {
    theResult = new G4FragmentVector();
    theResult->push_back(new G4Fragment(theFragment));
    return theResult;
  }

  G4double x = 1.0 + (etot - m0)/ekin;
  G4LorentzVector lv1;

  // scale and boost
  for (auto const & ptr : *theResult) {
    G4double m1 = ptr->GetGroundStateMass() + ptr->GetExcitationEnergy();
    G4double ek = std::max((ptr->GetMomentum().e() - m1)*x, 0.0);
    auto mom = ptr->GetMomentum().vect().unit();
    mom *= std::sqrt(ek * (ek + 2.0*m1));
    lv1.set(mom.x(), mom.y(), mom.z(), ek + m1);
    lv1.boost(bs);
    ptr->SetMomentum(lv1);
  }
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
  
  G4double T = std::max(Temperature, 0.0012*CLHEP::MeV);  
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

G4double G4StatMF::CalcEnergy(G4int A, G4int Z, const G4StatMFChannel* aChannel,
			      G4double T)
{
  G4double MassExcess0 = G4NucleiProperties::GetMassExcess(A,Z);
  G4double ChannelEnergy = aChannel->GetFragmentsEnergy(T);
  return -MassExcess0 + G4StatMFParameters::GetCoulomb() + ChannelEnergy;
}
