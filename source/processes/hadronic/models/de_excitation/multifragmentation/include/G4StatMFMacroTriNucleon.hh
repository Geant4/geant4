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
// $Id: G4StatMFMacroTriNucleon.hh,v 1.1 2003/08/26 18:47:26 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroTriNucleon_h
#define G4StatMFMacroTriNucleon_h 1

#include "G4VStatMFMacroCluster.hh"
#include "G4NucleiPropertiesTable.hh"

class G4StatMFMacroTriNucleon : public G4VStatMFMacroCluster {

public:

    // Default constructor
    G4StatMFMacroTriNucleon() : G4VStatMFMacroCluster(3) {};

    // Destructor
    ~G4StatMFMacroTriNucleon() {};
	

private:

    // Copy constructor
    G4StatMFMacroTriNucleon(const G4StatMFMacroTriNucleon & right);

    // operators
    G4StatMFMacroTriNucleon & operator=(const G4StatMFMacroTriNucleon & right);
    G4bool operator==(const G4StatMFMacroTriNucleon & right) const;
    G4bool operator!=(const G4StatMFMacroTriNucleon & right) const;

public:

    G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
				  const G4double nu, const G4double T);

    G4double CalcZARatio(const G4double ) {return theZARatio = 0.5;}

    G4double CalcEnergy(const G4double T);
	
    G4double CalcEntropy(const G4double T, const G4double FreeVol);
	
	
};

#endif
