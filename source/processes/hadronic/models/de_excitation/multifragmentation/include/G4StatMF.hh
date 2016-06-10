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
// $Id: G4StatMF.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4StatMF_h
#define G4StatMF_h 1

#include "globals.hh"
#include "G4VMultiFragmentation.hh"
#include "G4VStatMFEnsemble.hh"
#include "G4StatMFMicroCanonical.hh"
#include "G4StatMFMacroCanonical.hh"
#include "G4StatMFChannel.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

class G4StatMF : public G4VMultiFragmentation
{
public:
    // Default constructor
    G4StatMF();
    // Destructor
    ~G4StatMF();

private:
    // Copy constructor	
    G4StatMF(const G4StatMF & right);

    // Operators
    G4StatMF & operator=(const G4StatMF & right);
    G4bool operator==(const G4StatMF & right);
    G4bool operator!=(const G4StatMF & right);

public:

    G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);

private:

    // This finds temperature of breaking channel.
    G4bool FindTemperatureOfBreakingChannel(const G4Fragment & theFragment, 
					    const G4StatMFChannel * aChannel,
					    G4double & Temperature);

    // 
    G4double CalcEnergy(G4int A, G4int Z, 
			const G4StatMFChannel * aChannel,
			G4double T);


private:

    G4VStatMFEnsemble * _theEnsemble;



};

#endif












