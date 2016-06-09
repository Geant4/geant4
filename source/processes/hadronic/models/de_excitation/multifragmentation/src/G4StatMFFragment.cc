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
// $Id: G4StatMFFragment.cc,v 1.3 2003/11/04 11:30:35 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFFragment.hh"
#include "G4HadronicException.hh"


// Copy constructor
G4StatMFFragment::G4StatMFFragment(const G4StatMFFragment & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFFragment::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFFragment & G4StatMFFragment::
operator=(const G4StatMFFragment & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFFragment::operator==(const G4StatMFFragment & ) const
{
//	throw G4HadronicException(__FILE__, __LINE__, "G4StatMFFragment::operator== meant to not be accessable");
    return false;
}
 

G4bool G4StatMFFragment::operator!=(const G4StatMFFragment & ) const
{
//	throw G4HadronicException(__FILE__, __LINE__, "G4StatMFFragment::operator!= meant to not be accessable");
    return true;
}



G4double G4StatMFFragment::GetCoulombEnergy(void) const
{
    if (theZ <= 0.1) return 0.0;
    G4double Coulomb = (3./5.)*(elm_coupling*theZ*theZ)*
	pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.)/
	(G4StatMFParameters::Getr0()*pow(theA,1./3.));
						
    return Coulomb;
}


G4double G4StatMFFragment::GetEnergy(const G4double T) const
{
    if (theA < 1 || theZ < 0 || theZ > theA) {
	G4cerr << "G4StatMFFragment::GetEnergy: A = " << theA 
	       << ", Z = " << theZ << G4endl;
	throw G4HadronicException(__FILE__, __LINE__, 
	    "G4StatMFFragment::GetEnergy: Wrong values for A and Z!");
    }
    G4double BulkEnergy = G4NucleiProperties::GetMassExcess(static_cast<G4int>(theA),
							    static_cast<G4int>(theZ));
	
    if (theA < 4) return BulkEnergy - GetCoulombEnergy();
	
    G4double SurfaceEnergy;
    if (G4StatMFParameters::DBetaDT(T) == 0.0) SurfaceEnergy = 0.0;
    else SurfaceEnergy = (5./2.)*pow(theA,2.0/3.0)*T*T*
	     G4StatMFParameters::GetBeta0()/
	     (G4StatMFParameters::GetCriticalTemp()*
	      G4StatMFParameters::GetCriticalTemp());
					 
					 
    G4double ExchangeEnergy = theA*T*T/GetInvLevelDensity();
    if (theA != 4) ExchangeEnergy += SurfaceEnergy;		  
	
    return 	BulkEnergy + ExchangeEnergy - GetCoulombEnergy();		
	
}


G4double G4StatMFFragment::GetInvLevelDensity(void) const
{
    // Calculate Inverse Density Level
    // Epsilon0*(1 + 3 /(Af - 1))
    if (theA == 1) return 0.0;
    else return
	   G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(theA - 1.0));
}



G4Fragment * G4StatMFFragment::GetFragment(const G4double T)
{
    G4double U = CalcExcitationEnergy(T);
	
    G4double M = GetNuclearMass();

    G4LorentzVector FourMomentum(_momentum,sqrt(_momentum.mag2()+(M+U)*(M+U)));

    G4Fragment * theFragment = new G4Fragment(static_cast<G4int>(theA),static_cast<G4int>(theZ),FourMomentum);

    return theFragment;
}


G4double G4StatMFFragment::CalcExcitationEnergy(const G4double T)
{
    if (theA <= 3) return 0.0;
	
    G4double BulkEnergy = theA*T*T/GetInvLevelDensity();
	
    // if it is an alpha particle: done
    if (theA == 4) return BulkEnergy;
	
    // Term connected with surface energy
    G4double SurfaceEnergy = 0.0;
    if (abs(G4StatMFParameters::DBetaDT(T)) > 1.0e-20) 
// 		SurfaceEnergy = (5./2.)*pow(theA,2.0/3.0)*T*T*G4StatMFParameters::GetBeta0()/
// 			(G4StatMFParameters::GetCriticalTemp()*G4StatMFParameters::GetCriticalTemp());
	SurfaceEnergy = (5./2.)*pow(theA,2.0/3.0)*(G4StatMFParameters::Beta(T) - 
						   T*G4StatMFParameters::DBetaDT(T) - G4StatMFParameters::GetBeta0());
		
    return BulkEnergy + SurfaceEnergy;
}


