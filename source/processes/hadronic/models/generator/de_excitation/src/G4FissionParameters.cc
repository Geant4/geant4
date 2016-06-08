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
// $Id: G4FissionParameters.cc,v 1.5.2.1 2001/06/28 19:13:16 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4FissionParameters.hh"


const G4double G4FissionParameters::A1 = 134.0;
const G4double G4FissionParameters::A2 = 141.0;


G4FissionParameters::G4FissionParameters(const G4int A, const G4int Z, const G4double ExEnergy,
					 const G4double FissionBarrier)
{
    G4double U = ExEnergy; 
  
    As = A/2.0;

    if (A <= 235) Sigma2 = 5.6;  // MeV
    else Sigma2 = 5.6 + 0.096*(A-235); // MeV

    Sigma1 = 0.5*Sigma2; // MeV

    SigmaS = exp(0.00553*U/MeV + 2.1386); // MeV
    if (SigmaS > 20.0) SigmaS = 20.0;

    G4double FasymAsym = 2.0*exp(-((A2-As)*(A2-As))/(2.0*Sigma2*Sigma2)) + 
	exp(-((A1-As)*(A1-As))/(2.0*Sigma1*Sigma1));
  
    G4double FsymA1A2 = exp(-((As-(A1+A2)/2.0)*(As-(A1+A2)/2.0))/(2.0*SigmaS*SigmaS));
  

    G4double wa;
    G4double w1,w2;
    w = 0.0;
    if (Z >= 90) {         // Z >= 90
	if (U <= 16.25) wa = exp(0.5385*U/MeV-9.9564);  // U <= 16.25 MeV
	else wa = exp(0.09197*U/MeV-2.7003);            // U  > 16.25 MeV
    } else if (Z == 89) {  // Z == 89
	wa = exp(0.09197*U-1.0808);
    } else if (Z >= 82) {  //  82 <= Z <= 88
	G4double X = FissionBarrier - 7.5*MeV;
	if (X < 0.0) X = 0.0;
	wa = exp(0.09197*(U-X)/MeV-1.0808);
    } else {               // Z < 82
	w = 1001.0;
    }
  
    if (w == 0.0) {
	w1 = G4std::max(1.03*wa - FasymAsym, 0.0001);
	w2 = G4std::max(1.0 - FsymA1A2*wa,   0.0001);
    
	w = w1/w2;

	if (82 <= Z && Z < 89 && A < 227)  w *= exp(0.3*(227-A));
    }
  
}


G4FissionParameters::G4FissionParameters(const G4FissionParameters &right)
{
    G4Exception("G4FissionParameters::copy_constructor meant to not be accessable");
}


const G4FissionParameters & G4FissionParameters::operator=(const G4FissionParameters &right)
{
    G4Exception("G4FissionParameters::operator= meant to not be accessable");
    return *this;
}


G4bool G4FissionParameters::operator==(const G4FissionParameters &right) const
{
    return false;
}

G4bool G4FissionParameters::operator!=(const G4FissionParameters &right) const
{
    return true;
}


