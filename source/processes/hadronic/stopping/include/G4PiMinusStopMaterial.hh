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
// $Id: G4PiMinusStopMaterial.hh,v 1.11 2007-07-05 18:19:14 dennis Exp $
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
  virtual std::vector<G4ParticleDefinition*>* DefinitionVector();

  // 4-vectors of absorption products
  virtual std::vector<G4LorentzVector*>* P4Vector(const G4double binding, 
						   const G4double mass);

  // Number of final nucleons, out of generated absorption products
  virtual G4double FinalNucleons()=0;

protected:
   
  std::vector<G4ParticleDefinition* >* _definitions;
  std::vector<G4LorentzVector* >* _momenta; 
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
