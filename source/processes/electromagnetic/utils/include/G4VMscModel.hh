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
// $Id: G4VMscModel.hh,v 1.1 2003-07-21 12:59:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VMscModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 16.07.2003
//
// Modifications:
//
//

//
// Class Description:
//
// Abstract interface to models of multiple scattering

// -------------------------------------------------------------------
//

#ifndef G4VMscModel_h
#define G4VMscModel_h 1

#include "G4VEmModel.hh"

class G4PhysicsTable;

class G4VMscModel : public G4VEmModel
{

public:

  G4VMscModel(const G4String& val) : G4VEmModel(val) {};

  virtual ~G4VMscModel() {};

  virtual G4double GeomPathLength(G4PhysicsTable*,
                            const G4MaterialCutsCouple*,
		            const G4ParticleDefinition*,
		                  G4double&,
			          G4double,
			          G4double,
    			          G4double truePathLength) = 0;
  // G4double parameters: kinEnergy, lambda, range,
  // G4PhysicsTable: theLambdaTable

  virtual G4double TrueStepLength(G4double geomStepLength) = 0;

  virtual G4double SampleCosineTheta(G4double ) = 0;
  // G4double parameter trueStepLength

  virtual G4double SampleDisplacement() = 0;

private:

  //  hide assignment operator
  //  G4VMscModel & operator=(const  G4VMscModel &right);
  //  G4VMscModel(const  G4VMscModel&);
};

#endif

