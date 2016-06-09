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
// $Id: G4StatMFMacroNucleon.hh,v 1.8 2002/12/12 19:17:12 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroNucleon_h
#define G4StatMFMacroNucleon_h 1

#include "G4VStatMFMacroCluster.hh"


class G4StatMFMacroNucleon : public G4VStatMFMacroCluster {

public:

    // Default constructor
    G4StatMFMacroNucleon() : 
	G4VStatMFMacroCluster(1), _NeutronMeanMultiplicity(0.0),_ProtonMeanMultiplicity(0.0)
	 {};

    // Destructor
    ~G4StatMFMacroNucleon() {};
	

private:

    // Copy constructor
    G4StatMFMacroNucleon(const G4StatMFMacroNucleon & right);

    // operators
    G4StatMFMacroNucleon & operator=(const G4StatMFMacroNucleon & right);
    G4bool operator==(const G4StatMFMacroNucleon & right) const;
    G4bool operator!=(const G4StatMFMacroNucleon & right) const;

public:

    G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
				  const G4double nu, const G4double T);
	
    G4double CalcZARatio(const G4double nu) 
	{ if (_ProtonMeanMultiplicity+_NeutronMeanMultiplicity > 0.0)
	    return theZARatio = _ProtonMeanMultiplicity/
		(_ProtonMeanMultiplicity+_NeutronMeanMultiplicity);
	else return 0.0; }
						
						
    G4double CalcEnergy(const G4double T);
	
    G4double CalcEntropy(const G4double T, const G4double FreeVol);

	
private:

    G4double _NeutronMeanMultiplicity;
    G4double _ProtonMeanMultiplicity;

};

#endif
