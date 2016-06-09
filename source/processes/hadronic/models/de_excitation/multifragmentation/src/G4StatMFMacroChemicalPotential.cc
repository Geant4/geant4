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
// $Id: G4StatMFMacroChemicalPotential.cc,v 1.4 2005/06/04 13:27:48 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFMacroChemicalPotential.hh"

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
    G4double CP = ((3./5.)*elm_coupling/G4StatMFParameters::Getr0())*
	(1.0-1.0/std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0));

    // Initial value for _ChemPotentialNu	
    _ChemPotentialNu = (theZ/theA)*(8.0*G4StatMFParameters::GetGamma0()+2.0*CP*std::pow(theA,2./3.)) -
	4.0*G4StatMFParameters::GetGamma0();
		

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
	    } while (fChemPb < 0.0);
	} else {
	    do {
		ChemPb += 1.5*std::abs(ChemPb-ChemPa);
		fChemPb = this->operator()(ChemPb);
	    } while (fChemPb > 0.0);
	}
    }

    G4Solver<G4StatMFMacroChemicalPotential> * theSolver =
      new G4Solver<G4StatMFMacroChemicalPotential>(100,1.e-4);
    theSolver->SetIntervalLimits(ChemPa,ChemPb);
    //    if (!theSolver->Crenshaw(*this)) 
    if (!theSolver->Brent(*this)) 
      throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMacroChemicalPotential::CalcChemicalPotentialNu: I couldn't find the root.");
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
      MeanZ += static_cast<G4double>(n++) *
	(*i)->GetZARatio() *
	(*i)->GetMeanMultiplicity(); 
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
