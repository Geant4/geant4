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
// $Id: G4PiMinusStopMaterial.hh,v 1.7 2001-10-04 20:00:40 hpw Exp $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PiMinusStopMaterial.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 18 May 1998
//
//      Modifications: 
//      13 Sep 1998 - MGP   Modified P4Vector
//      
// -------------------------------------------------------------------

#ifndef G4PIMINUSSTOPMATERIAL_HH
#define G4PIMINUSSTOPMATERIAL_HH

#include "globals.hh"
#include "G4LorentzVector.hh"
//#include "G4String.hh"
#include "G4DistributionGenerator.hh"
#include "G4ParticleDefinition.hh"


class G4PiMinusStopMaterial 
{  

private:

  // Hide assignment operator as private 
  G4PiMinusStopMaterial& operator=(const G4PiMinusStopMaterial &right);

  // Copy constructor
  G4PiMinusStopMaterial(const G4PiMinusStopMaterial& );

public:

  // Constructor
  G4PiMinusStopMaterial();

  // Destructor
  virtual ~G4PiMinusStopMaterial();

  // Definitions of absorption products
  virtual G4std::vector<G4ParticleDefinition*>* DefinitionVector();

  // 4-vectors of absorption products
  virtual G4std::vector<G4LorentzVector*>* P4Vector(const G4double binding, 
						   const G4double mass);

  // Number of final nucleons, out of generated absorption products
  virtual G4double FinalNucleons()=0;

protected:
   
  G4std::vector<G4ParticleDefinition* >* _definitions;
  G4std::vector<G4LorentzVector* >* _momenta; 
  G4DistributionGenerator* _distributionE;
  G4DistributionGenerator* _distributionAngle;
  G4double theR; 

  G4double GenerateAngle(G4double range);
  G4LorentzVector MakeP4(G4double p, G4double theta, 
			 G4double phi, G4double e);
  G4double RecoilEnergy(const G4double mass);

private:

};
 

#endif
