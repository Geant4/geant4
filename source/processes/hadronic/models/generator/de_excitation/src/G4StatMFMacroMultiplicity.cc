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
// $Id: G4StatMFMacroMultiplicity.cc,v 1.11 2003/06/30 07:43:58 gcosmo Exp $
// GEANT4 tag $Name: ghad-deex-V05-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFMacroMultiplicity.hh"

// operators definitions
G4StatMFMacroMultiplicity & 
G4StatMFMacroMultiplicity::operator=(const G4StatMFMacroMultiplicity & ) 
{
    G4Exception("G4StatMFMacroMultiplicity::operator= meant to not be accessable");
    return *this;
}

G4bool G4StatMFMacroMultiplicity::operator==(const G4StatMFMacroMultiplicity & ) const 
{
    G4Exception("G4StatMFMacroMultiplicity::operator== meant to not be accessable");
    return false;
}


G4bool G4StatMFMacroMultiplicity::operator!=(const G4StatMFMacroMultiplicity & ) const 
{
    G4Exception("G4StatMFMacroMultiplicity::operator!= meant to not be accessable");
    return true;
}




G4double G4StatMFMacroMultiplicity::CalcChemicalPotentialMu(void) 
    //	Calculate Chemical potential \mu
    // For that is necesary to calculate mean multiplicities
{
    G4double CP = ((3./5.)*elm_coupling/G4StatMFParameters::Getr0())*
	(1.0-1.0/pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1.0/3.0));

    // starting value for chemical potential \mu
    // it is the derivative of F(T,V)-\nu*Z w.r.t. Af in Af=5
    G4double ZA5 = _theClusters->operator[](4)->GetZARatio();
    G4double ILD5 = _theClusters->operator[](4)->GetInvLevelDensity();
    _ChemPotentialMu = -G4StatMFParameters::GetE0()-
	_MeanTemperature*_MeanTemperature/ILD5 -
	_ChemPotentialNu*ZA5 + 
	G4StatMFParameters::GetGamma0()*(1.0-2.0*ZA5)*(1.0-2.0*ZA5) +
	(2.0/3.0)*G4StatMFParameters::Beta(_MeanTemperature)/pow(5.,1./3.) +
	(5.0/3.0)*CP*ZA5*ZA5*pow(5.,2./3.) -
	1.5*_MeanTemperature/5.0;
		


    G4double ChemPa = _ChemPotentialMu;
    if (ChemPa/_MeanTemperature > 10.0) ChemPa = 10.0*_MeanTemperature;
    G4double ChemPb = ChemPa - 0.5*abs(ChemPa);
    
    
    G4double fChemPa = this->operator()(ChemPa); 
    G4double fChemPb = this->operator()(ChemPb); 
    
    // bracketing the solution
    G4int iterations = 0;
    while (fChemPa*fChemPb > 0.0 && iterations < 10) 
    {
	if (abs(fChemPa) <= abs(fChemPb)) 
	{
	    ChemPa += 0.6*(ChemPa-ChemPb);
	    fChemPa = this->operator()(ChemPa);
	} 
	else 
	{
	    ChemPb += 0.6*(ChemPb-ChemPa);
	    fChemPb = this->operator()(ChemPb);
	}
    }
    if (fChemPa*fChemPb > 0.0) 
    {
	G4Exception("G4StatMFMacroMultiplicity::CalcChemicalPotentialMu: I couldn't bracket the root.");
    }
	
	
    G4Solver<G4StatMFMacroMultiplicity> * theSolver = new G4Solver<G4StatMFMacroMultiplicity>(100,1.e-4);
    theSolver->SetIntervalLimits(ChemPa,ChemPb);
    //    if (!theSolver->Crenshaw(*this)) 
    if (!theSolver->Brent(*this)) 
    {
	G4Exception("G4StatMFMacroMultiplicity::CalcChemicalPotentialMu: I couldn't find the root.");
    }
    _ChemPotentialMu = theSolver->GetRoot();
    delete theSolver;
    return _ChemPotentialMu;
}



G4double G4StatMFMacroMultiplicity::CalcMeanA(const G4double mu)
{
  G4double r03 = G4StatMFParameters::Getr0(); r03 *= r03*r03;
  G4double V0 = (4.0/3.0)*pi*theA*r03;

  G4double MeanA = 0.0;
	
  _MeanMultiplicity = 0.0;
	
 
  G4int n = 1;
 for (std::vector<G4VStatMFMacroCluster*>::iterator i = _theClusters->begin(); 
      i != _theClusters->end(); ++i) 
   {
     G4double multip = (*i)->CalcMeanMultiplicity(V0*_Kappa,mu,_ChemPotentialNu,_MeanTemperature);
     MeanA += multip*static_cast<G4double>(n++);
     _MeanMultiplicity += multip;
   }

  return MeanA;
}
