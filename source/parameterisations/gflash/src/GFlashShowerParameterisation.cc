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
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ------- GFlashShowerParameterisation -------
//
// Authors: Joanna Weng - 11.2005
// ------------------------------------------------------------

#include "GFlashShowerParameterisation.hh"
#include <cmath>
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Material.hh"
#include "Gamma.hh" // @@@@
#include "G4MaterialTable.hh"

GFlashShowerParameterisation::GFlashShowerParameterisation()
{  
	
	
}

G4double GFlashShowerParameterisation::GetEffZ(const G4Material * mat  )
{
	// Returns Z or effective Z=sum(pi*Zi) (if compound/mixture)
	// of given material.
	G4double z = 0.;
	G4int nofElements = mat->GetNumberOfElements();
	if (nofElements > 1) 
	{
		for (G4int i=0; i<nofElements; i++) {
			G4double zOfElement = mat->GetElement(i)->GetZ();
			G4double massFraction = mat->GetFractionVector()[i];
			// cout << mat->GetElement(i)->GetName() <<" Z= "<<zOfElement << " , Fraction= "<<massFraction <<endl;
			z += zOfElement*massFraction;
		}
	}
	else { 
		z = mat->GetZ(); 
	}  
	return z;
}


G4double GFlashShowerParameterisation::GetEffA  (const G4Material * mat  )
{
	// Returns A or effective A=sum(pi*Ai) (if compound/mixture)
	// of given material.
	G4double a = 0.;
	G4int nofElements = mat->GetNumberOfElements();
	if (nofElements > 1) {
		for (G4int i=0; i<nofElements; i++) {
			G4double aOfElement = mat->GetElement(i)->GetA()/(g/mole);
			G4double massFraction = mat->GetFractionVector()[i];     
			a += aOfElement*massFraction;
		}
	}
	else { 
		a = mat->GetA()/(g/mole);
	}
	return a;
}

void GFlashShowerParameterisation::PrintMaterial(const G4Material * mat)
{
	G4cout<<"/********************************************/ " << G4endl;
	G4cout<<"  - GFlashShowerParameterisation::Material -  " << G4endl;
	G4cout<<"        Material : " << mat->GetName()  << G4endl;
	G4cout<<"   Z = "<< Z  << G4endl;
	G4cout<<"   A = "<< A  << G4endl;
	G4cout<<"   X0 = "<<X0/cm <<" cm" << G4endl;
	G4cout<<"    Rm= "<<Rm/cm <<" cm" << G4endl;
	G4cout<<"   Ec = "<<Ec/MeV << " MeV"<< G4endl;
	G4cout<<"/********************************************/ " << G4endl; 
}

GFlashShowerParameterisation::~GFlashShowerParameterisation()
{}


G4double GFlashShowerParameterisation::GeneratePhi()
{
	G4double Phi = twopi*G4UniformRand() ;
	return Phi;
}

G4double GFlashShowerParameterisation::gam(G4double x, G4double a) const 
{
	static MyGamma theG;
	return  theG.Gamma(a, x); 
}
