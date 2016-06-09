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
// $Id: G4StatMFFragment.cc,v 1.7 2008-07-25 11:20:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
	std::pow(1.0+G4StatMFParameters::GetKappaCoulomb(),1./3.)/
	(G4StatMFParameters::Getr0()*std::pow(theA,1./3.));
						
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
    else SurfaceEnergy = (5./2.)*std::pow(theA,2.0/3.0)*T*T*
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

    G4LorentzVector FourMomentum(_momentum,std::sqrt(_momentum.mag2()+(M+U)*(M+U)));

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
    if (std::abs(G4StatMFParameters::DBetaDT(T)) > 1.0e-20) 
// 		SurfaceEnergy = (5./2.)*std::pow(theA,2.0/3.0)*T*T*G4StatMFParameters::GetBeta0()/
// 			(G4StatMFParameters::GetCriticalTemp()*G4StatMFParameters::GetCriticalTemp());
	SurfaceEnergy = (5./2.)*std::pow(theA,2.0/3.0)*(G4StatMFParameters::Beta(T) - 
						   T*G4StatMFParameters::DBetaDT(T) - G4StatMFParameters::GetBeta0());
		
    return BulkEnergy + SurfaceEnergy;
}


