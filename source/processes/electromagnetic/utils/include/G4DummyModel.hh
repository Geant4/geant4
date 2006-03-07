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
// $Id: G4DummyModel.hh,v 1.1 2006-03-07 16:57:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4DummyModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2006
//
// Modifications:
//
// Class Description:
//
// EM model doing nothing

// -------------------------------------------------------------------
//

#ifndef G4DummyModel_h
#define G4DummyModel_h 1

#include "globals.hh"
#include "G4VEmModel.hh"

class G4DummyModel :  public G4VEmModel
{

public:

  G4DummyModel(const G4String& nam = "DummyModel");

  virtual ~G4DummyModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax);

private: 

  //  hide assignment operator
  G4DummyModel & operator=(const  G4DummyModel &right);
  G4DummyModel(const  G4DummyModel&);
};

#endif

