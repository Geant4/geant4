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
// $Id: G4StatMFMacroMultiNucleon.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroMultiNucleon_h
#define G4StatMFMacroMultiNucleon_h 1

#include "G4VStatMFMacroCluster.hh"


class G4StatMFMacroMultiNucleon : public G4VStatMFMacroCluster {

public:

    // Constructor
    G4StatMFMacroMultiNucleon(const G4int Size) : G4VStatMFMacroCluster(Size) {};

    // Destructor
    ~G4StatMFMacroMultiNucleon() {};
	

private:

    // Default constructor
    G4StatMFMacroMultiNucleon();
	
    // Copy constructor
    G4StatMFMacroMultiNucleon(const G4StatMFMacroMultiNucleon & right);

    // operators
    G4StatMFMacroMultiNucleon & operator=(const G4StatMFMacroMultiNucleon & right);
    G4bool operator==(const G4StatMFMacroMultiNucleon & right) const;
    G4bool operator!=(const G4StatMFMacroMultiNucleon & right) const;

public:

    G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
				  const G4double nu, const G4double T);
								
    G4double CalcZARatio(const G4double nu); 
	
    G4double CalcEnergy(const G4double T);

    G4double CalcEntropy(const G4double T, const G4double FreeVol);
	
};

#endif
