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
// $Id: G4PenelopeBremsstrahlungFS.hh 75573 2013-11-04 11:48:15Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 20 Oct 2010   L. Pandola   1st implementation. 
// 02 May 2011   L. Pandola   Remove dependency on CLHEP::HepMatrix
// 24 May 2011   L. Pandola   Renamed (make default Penelope)
// 03 Oct 2013   L. Pandola   Migration to MT
// 07 Oct 2013   L. Pandola   Add verbosity and ismaster flag for the 
//                             master-only methods
// 30 Oct 2013   L. Pandola   Add a G4Cache of temporary vectors as 
//                             private member. 
//
// -------------------------------------------------------------------
//
// Class description:
// Helper class for the calculation of final state (energy and angular 
// distribution) for bremsstrahlung, Penelope Model, version 2008
// -------------------------------------------------------------------

#ifndef G4PENELOPEBREMSSTRAHLUNGFS_HH
#define G4PENELOPEBREMSSTRAHLUNGFS_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4Cache.hh"
#include <map>

class G4PhysicsFreeVector;
class G4PhysicsLogVector;
class G4Material;
class G4PhysicsTable;

class G4PenelopeBremsstrahlungFS 
{
public:
  //! Only master models are supposed to create instances
  G4PenelopeBremsstrahlungFS(G4int verbosity=0);
  ~G4PenelopeBremsstrahlungFS();
 
  //!
  //! Master and workers (do not touch tables)
  //! All of them are const
  //!
  G4double GetEffectiveZSquared(const G4Material* mat) const;
  size_t GetNBinsX() const {return nBinsX;};
  G4double GetMomentumIntegral(G4double* y,			   
			       G4double up,G4int momOrder) const;
  const G4PhysicsTable* GetScaledXSTable(const G4Material*,
					 const G4double cut) const;
  G4double SampleGammaEnergy(G4double energy,
			     const G4Material*, const G4double cut) const;

  //! Reserved for the master model: they build and handle tables
  void ClearTables(G4bool isMaster=true);
  void BuildScaledXSTable(const G4Material* material,G4double cut,
			  G4bool isMaster);

  void SetVerbosity(G4int ver){fVerbosity=ver;};
  G4int GetVerbosity(){return fVerbosity;};

private:
  //assignment operator
  G4PenelopeBremsstrahlungFS & operator=(const G4PenelopeBremsstrahlungFS &right);
  //copy constructor
  G4PenelopeBremsstrahlungFS(const G4PenelopeBremsstrahlungFS&);
  
  //Differential cross section tables
  //Table contains G4PhysicsVectors of log(XS(E,x)) vs. log(E) 
  //for a grid of 32 values in x (= reduced photon energy)
  std::map< std::pair<const G4Material*,G4double> , 
	    G4PhysicsTable*> *theReducedXSTable; 

  std::map<const G4Material*,G4double> *theEffectiveZSq;
  

  //Element data table
  static const size_t nBinsE = 57;
  static const size_t nBinsX = 32;
  //x and E grids used in the data tables
  G4double theXGrid[nBinsX]; 
  G4double theEGrid[nBinsE];
  void ReadDataFile(G4int Z);  

  //Map of element data vs. Z. 
  //This is conceptually an array [57][33] with 57 energy 
  //points and 32 points in x. The 33-th column gives the total XS vs. E.
  //It is implemented as a one-dimensional array of dimension
  //57*33=1881 elements. data[e][x] --> theElementData[e*(nBinsX+1)+x]
  std::map<G4int,G4DataVector*> *theElementData;

  //Tables for energy sampling
  void InitializeEnergySampling(const G4Material*,G4double cut);
  
  //Table contains G4PhysicsVectors of integral_XS(E,x) vs. x for a grid of 
  //57 values in energy
  std::map< std::pair<const G4Material*,G4double> , 
	    G4PhysicsTable*> *theSamplingTable;
  std::map< std::pair<const G4Material*,G4double> , 
	    G4PhysicsFreeVector* > *thePBcut;
 
  //temporary vector. Used as member variable to avoid to book/release the 
  //memory on the fly. This vector is over-written at every call of 
  //SampleGammaEnergy(). It is thread-local (each thread has its own) 
  //and managed by G4Cache
  G4Cache<G4PhysicsFreeVector*> fCache;

  G4int fVerbosity;

};

#endif

