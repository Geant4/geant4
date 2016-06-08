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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMFMacroTemperature.cc,v 1.8 2001/08/01 17:05:34 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFMacroTemperature.hh"

// operators definitions
G4StatMFMacroTemperature & 
G4StatMFMacroTemperature::operator=(const G4StatMFMacroTemperature & right) 
{
    G4Exception("G4StatMFMacroTemperature::operator= meant to not be accessable");
    return *this;
}

G4bool G4StatMFMacroTemperature::operator==(const G4StatMFMacroTemperature & right) const 
{
    G4Exception("G4StatMFMacroTemperature::operator== meant to not be accessable");
    return false;
}


G4bool G4StatMFMacroTemperature::operator!=(const G4StatMFMacroTemperature & right) const 
{
    G4Exception("G4StatMFMacroTemperature::operator!= meant to not be accessable");
    return true;
}




G4double G4StatMFMacroTemperature::CalcTemperature(void) 
    //	Calculate Chemical potential \nu
{
    // Temperature
    G4double Ta = 0.00012; 
    G4double Tb = G4std::max(sqrt(_ExEnergy/(theA*0.12)),0.01*MeV);
    
    G4double fTa = this->operator()(Ta); 
    G4double fTb = this->operator()(Tb); 

    // Bracketing the solution
    // T should be greater than 0.
    // The interval is [Ta,Tb]
    // We start with a low value for Ta = 0.0012 K
    // it should be enough to have fTa > 0 If it isn't 
    // the case, we decrease Ta. But carefully, because 
    // fTa growes very fast when Ta is near 0 and we could have
    // an overflow.

    G4int iterations = 0;  
    while (fTa < 0.0 && iterations++ < 10) {
	Ta -= 0.5*Ta;
	fTa = this->operator()(Ta);
    }
    // Usually, fTb will be less than 0, but if it is not the case: 
    iterations = 0;  
    while (fTa*fTb > 0.0 && iterations++ < 10) {
	Tb += 1.5*abs(Tb-Ta);
	fTb = this->operator()(Tb);
    }
	
    if (fTa*fTb > 0.0) {
	G4Exception("G4StatMFMacroTemperature::CalcTemperature: I couldn't bracket	the solution.");
    }

    G4Solver<G4StatMFMacroTemperature> * theSolver = new G4Solver<G4StatMFMacroTemperature>(100,1.e-4);
    theSolver->SetIntervalLimits(Ta,Tb);
    if (!theSolver->Brent(*this)) 
	G4Exception("G4StatMFMacroTemperature::CalcTemperature: I couldn't find the root.");
    _MeanTemperature = theSolver->GetRoot();
    delete theSolver;
    return _MeanTemperature;
}



G4double G4StatMFMacroTemperature::FragsExcitEnergy(const G4double T)
    // Calculates excitation energy per nucleon and summed fragment multiplicity and entropy
{

    // Model Parameters
    G4double R0 = G4StatMFParameters::Getr0()*pow(theA,1./3.);
    G4double R = R0*pow(1.0+G4StatMFParameters::GetKappaCoulomb(), 1./3.);
    G4double FreeVol = _Kappa*(4.*pi/3.)*R0*R0*R0; 
 
 
    // Calculate Chemical potentials
    CalcChemicalPotentialNu(T);


    // Average total fragment energy
    G4double AverageEnergy = 0.0;
    G4int i;
    for (i = 0; i < theA; i++) AverageEnergy += 
				   _theClusters->operator[](i)->GetMeanMultiplicity()*
				   _theClusters->operator[](i)->CalcEnergy(T);


    // Add Coulomb energy			
    AverageEnergy += (3./5.)*elm_coupling*theZ*theZ/R;		

    // Calculate mean entropy
    _MeanEntropy = 0.0;
    for (i = 0; i < theA; i++) _MeanEntropy +=
				   _theClusters->operator[](i)->CalcEntropy(T,FreeVol);	

    // Excitation energy per nucleon
    G4double FragsExcitEnergy = AverageEnergy - _FreeInternalE0;

    return FragsExcitEnergy;

}


void G4StatMFMacroTemperature::CalcChemicalPotentialNu(const G4double T)
    // Calculates the chemical potential \nu 

{
    G4StatMFMacroChemicalPotential * theChemPot = new
	G4StatMFMacroChemicalPotential(theA,theZ,_Kappa,T,_theClusters);


    _ChemPotentialNu = theChemPot->CalcChemicalPotentialNu();
    _ChemPotentialMu = theChemPot->GetChemicalPotentialMu();
    _MeanMultiplicity = theChemPot->GetMeanMultiplicity();	
	
    delete theChemPot;
		    
    return;

}


