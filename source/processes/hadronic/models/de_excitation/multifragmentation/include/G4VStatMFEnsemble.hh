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
// $Id: G4VStatMFEnsemble.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4VStatMFEnsemble_h
#define G4VStatMFEnsemble_h 1

#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

class G4VStatMFEnsemble {

public:
    // Default Constructor
    G4VStatMFEnsemble() :
	__FreeInternalE0(0.0),
	__MeanTemperature(0.0),
	__MeanEntropy(0.0),
	__MeanMultiplicity(0.0)
	{};


    // Destructor
    virtual ~G4VStatMFEnsemble() {};


private:

    // Copy constructor
    G4VStatMFEnsemble(const G4VStatMFEnsemble & right);

    // operators
    G4VStatMFEnsemble & operator=(const G4VStatMFEnsemble & right);
    G4bool operator==(const G4VStatMFEnsemble & right) const;
    G4bool operator!=(const G4VStatMFEnsemble & right) const;
	
public:

    virtual G4StatMFChannel * ChooseAandZ(const G4Fragment & aFragment) = 0;
		
    G4double GetMeanMultiplicity(void) const {return __MeanMultiplicity;}
	
    G4double GetMeanTemperature(void) const {return __MeanTemperature;}

protected:

    // Free internal energy at temperature T = 0
    G4double __FreeInternalE0;
	
    // Mean temperature 
    G4double __MeanTemperature;
	
    // Mean Entropy 
    G4double __MeanEntropy;
	
    // Mean Multiplicity
    G4double __MeanMultiplicity;
};

#endif
