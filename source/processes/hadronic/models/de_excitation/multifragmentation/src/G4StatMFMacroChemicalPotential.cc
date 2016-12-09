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
// $Id: G4StatMFMacroChemicalPotential.cc 100379 2016-10-19 15:05:35Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMacroChemicalPotential.hh"
#include "G4PhysicalConstants.hh"
#include "G4Pow.hh"

// operators definitions
G4StatMFMacroChemicalPotential & 
G4StatMFMacroChemicalPotential::operator=(const G4StatMFMacroChemicalPotential & ) 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroChemicalPotential::operator= meant to not be accessable");
    return *this;
}

G4bool G4StatMFMacroChemicalPotential::operator==(const G4StatMFMacroChemicalPotential & ) const 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroChemicalPotential::operator== meant to not be accessable");
    return false;
}


G4bool G4StatMFMacroChemicalPotential::operator!=(const G4StatMFMacroChemicalPotential & ) const 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroChemicalPotential::operator!= meant to not be accessable");
    return true;
}

G4double G4StatMFMacroChemicalPotential::CalcChemicalPotentialNu(void) 
//	Calculate Chemical potential \nu
{
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double CP = G4StatMFParameters::GetCoulomb();

  // Initial value for _ChemPotentialNu	
  _ChemPotentialNu = (theZ/theA)*(8.0*G4StatMFParameters::GetGamma0()
				  +2.0*CP*g4calc->Z23(theA))
    - 4.0*G4StatMFParameters::GetGamma0();
		
  G4double ChemPa = _ChemPotentialNu;
  G4double ChemPb = 0.5*_ChemPotentialNu;
    
  G4double fChemPa = this->operator()(ChemPa); 
  G4double fChemPb = this->operator()(ChemPb); 

  if (fChemPa*fChemPb > 0.0) {    
    // bracketing the solution
    if (fChemPa < 0.0) {
      do {
	ChemPb -= 1.5*std::abs(ChemPb-ChemPa);
	fChemPb = this->operator()(ChemPb);   
	// Loop checking, 05-Aug-2015, Vladimir Ivanchenko
      } while (fChemPb < 0.0);
    } else {
      do {
	ChemPb += 1.5*std::abs(ChemPb-ChemPa);
	fChemPb = this->operator()(ChemPb);
	// Loop checking, 05-Aug-2015, Vladimir Ivanchenko
      } while (fChemPb > 0.0);
    }
  }

  G4Solver<G4StatMFMacroChemicalPotential> * theSolver =
    new G4Solver<G4StatMFMacroChemicalPotential>(100,1.e-4);
  theSolver->SetIntervalLimits(ChemPa,ChemPb);
  //    if (!theSolver->Crenshaw(*this)) 
  if (!theSolver->Brent(*this)){
    G4cout <<"G4StatMFMacroChemicalPotential:"<<" ChemPa="<<ChemPa
	   <<" ChemPb="<<ChemPb<< G4endl;
    G4cout <<"G4StatMFMacroChemicalPotential:"<<" fChemPa="<<fChemPa
	   <<" fChemPb="<<fChemPb<< G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroChemicalPotential::CalcChemicalPotentialNu: I couldn't find the root.");
  }
  _ChemPotentialNu = theSolver->GetRoot();
  delete theSolver;
  return _ChemPotentialNu;
}



G4double G4StatMFMacroChemicalPotential::CalcMeanZ(const G4double nu)
{
  std::vector<G4VStatMFMacroCluster*>::iterator i;
  for (i= _theClusters->begin()+1; i != _theClusters->end(); ++i) 
    { 
      (*i)->CalcZARatio(nu);
    }
  CalcChemicalPotentialMu(nu);
  // This is important, the Z over A ratio for proton and neutron depends on the 
  // chemical potential Mu, while for the first guess for Chemical potential mu 
  // some values of Z over A ratio. This is the reason for that.
  (*_theClusters->begin())->CalcZARatio(nu);
  
  G4double MeanZ = 0.0;
  G4int n = 1;
  for (i = _theClusters->begin(); i != _theClusters->end(); ++i)
    {
      MeanZ += (n++) * (*i)->GetZARatio() * (*i)->GetMeanMultiplicity(); 
    }
  return MeanZ;
}

void G4StatMFMacroChemicalPotential::CalcChemicalPotentialMu(const G4double nu)
//	Calculate Chemical potential \mu
// For that is necesary to calculate mean multiplicities
{
  G4StatMFMacroMultiplicity * theMultip = new
    G4StatMFMacroMultiplicity(theA,_Kappa,_MeanTemperature,nu,_theClusters);
  
  _ChemPotentialMu = theMultip->CalcChemicalPotentialMu();
  _MeanMultiplicity = theMultip->GetMeanMultiplicity();
  
  delete theMultip;
  
  return; 
}
