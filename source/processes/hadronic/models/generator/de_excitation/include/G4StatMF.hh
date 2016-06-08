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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMF.hh,v 1.7 2002/01/15 12:18:42 vlara Exp $
// GEANT4 tag $Name: geant4-04-01 $
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


//#define pctest

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
					    const G4StatMFChannel * aChannel,G4double & Temperature);

    // 
    G4double CalcEnergy(const G4double A, const G4double Z, 
			const G4StatMFChannel * aChannel,
			const G4double T);


private:

    G4VStatMFEnsemble * _theEnsemble;



};

#endif












