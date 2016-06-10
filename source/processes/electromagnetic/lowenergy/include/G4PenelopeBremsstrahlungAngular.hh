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
// -------------------------------------------------------------------
// $Id: G4PenelopeBremsstrahlungAngular.hh 74452 2013-10-07 15:08:00Z gcosmo $
//
// Author: L.Pandola
//
// History:
// -----------
// 23 Nov 2010  L. Pandola       1st implementation
// 24 May 2011  L. Pandola       Renamed (make default Penelope)
// 13 Mar 2012  L. Pandola       Made a derived class of G4VEmAngularDistribution
//                               and update the interface accordingly
// 18 Jul 2012  L. Pandola       Migrate to the new interface of G4VEmAngularDistribution
// 03 Oct 2013  L. Pandola       Migration to MT
//
// Class description:
// Calculation of angular distribution for Penelope Bremsstrahlung
// version 2008
// --------------------------------------------------------------


#ifndef G4PENELOPEBREMSSTRAHLUNGANGULAR_HH
#define G4PENELOPEBREMSSTRAHLUNGANGULAR_HH 1
#include "globals.hh"
#include <map>
#include "G4VEmAngularDistribution.hh"

class G4PhysicsTable;
class G4Material;

class G4PenelopeBremsstrahlungAngular : public G4VEmAngularDistribution
{ 

public:
  G4PenelopeBremsstrahlungAngular(); 
  ~G4PenelopeBremsstrahlungAngular();

  
  //! Old interface, backwards compatibility. Will not work in this case
  //! it will produce a G4Exception().
  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  //! Samples the direction of the outgoing photon (in global coordinates). 
  G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
				 G4double out_energy,
				 G4int Z,
				 const G4Material* mat = 0);
  
  //! Set/Get Verbosity level
  void SetVerbosityLevel(G4int vl){verbosityLevel = vl;};
  G4int GetVerbosityLevel(){return verbosityLevel;};

  //! Reserved for Master Model
  //! The Initialize() method forces the cleaning of tables
  void Initialize();
  //! Reserved for Master Model
  void PrepareTables(const G4Material* material,
		     G4bool isMaster);


private:
  
  void ClearTables();
 
  G4double CalculateEffectiveZ(const G4Material* material);
  
  std::map<const G4Material*,G4double> *theEffectiveZSq;

  //Tables containing the Lorentz sampling coefficients 
  //The key is the effective Z of the material
  std::map<G4double,G4PhysicsTable*> *theLorentzTables1;
  std::map<G4double,G4PhysicsTable*> *theLorentzTables2;

  void ReadDataFile();
  G4bool dataRead;
  
  static const G4int NumberofZPoints=6;
  static const G4int NumberofEPoints=6;
  static const G4int NumberofKPoints=4;

  G4double QQ1[NumberofZPoints][NumberofEPoints][NumberofKPoints];
  G4double QQ2[NumberofZPoints][NumberofEPoints][NumberofKPoints];

  G4int verbosityLevel;

};


  
#endif
