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
// $Id: G4VTRModel.hh,v 1.2 2005-12-19 15:08:41 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VTRModel  -- header file
//
// The model of transition radiation
//
// History:
//
// 04.10.05, V.Grichine move from pure virtual and new class name 
// 29.02.04, V.Ivanchenko created

#ifndef G4VTRModel_h
#define G4VTRModel_h


#include "globals.hh"
#include <vector>

class G4Material;
class G4Track;
class G4VParticleChange;


class G4VTRModel
{
public:

// Constructors

  G4VTRModel( const G4String& modelName) {fName = modelName;};

// Destructor

  virtual ~G4VTRModel() {};

  const G4String& GetName() const {return fName;};

  virtual void GenerateSecondaries(G4VParticleChange& pChange,
                                   std::vector<const G4Material*>& materials,
                                   std::vector<G4double>& steps,
                                   std::vector<G4ThreeVector>& normals,
		                   G4ThreeVector& startingPosition,
		             const G4Track& track);

  virtual void PrintInfo() { return; };

// private :

  // hide assignment operator

  G4VTRModel & operator=(const G4VTRModel &right);
  G4VTRModel(const G4VTRModel&);

protected:


  G4String  fName;
 

};

#endif   // G4VTRModel_h
