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
/*
 * File:   G4WendtFissionFragmentGenerator.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 21, 2013, 13:58 MST
 */

#ifndef G4WendtFissionFragmentGenerator_h
#define G4WendtFissionFragmentGenerator_h 1

#include <map>

#include "G4HadFinalState.hh"

#include "G4FissionFragmentGenerator.hh"

class G4WendtFissionFragmentGenerator
{
public:
	G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile, G4int Z, G4int A);
	static G4WendtFissionFragmentGenerator* GetInstance() {
           if ( instance == NULL) instance = new G4WendtFissionFragmentGenerator();
           return instance;
        }

	void InitializeANucleus(const G4int A, const G4int Z, const G4int M, const G4String& dataDirectory);
	~G4WendtFissionFragmentGenerator();

private:
	// SINGLETON!!!
	G4WendtFissionFragmentGenerator();
	G4WendtFissionFragmentGenerator(G4WendtFissionFragmentGenerator const&);
	void operator=(G4WendtFissionFragmentGenerator const&);

        static G4ThreadLocal G4WendtFissionFragmentGenerator* instance;
	// SINGLETON!!!

	/** A map of all the fission isotopes loaded at initialization */
    std::map< const G4int, G4FissionFragmentGenerator* > fissionIsotopes;
    G4ParticleHPNames fileNames;

    G4int Verbosity_;
};
#endif
