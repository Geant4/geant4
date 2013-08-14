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
// $Id: G4MultiBodyMomentumDist.hh 71719 2013-06-21 00:01:54Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: Singleton class to evaluate multi-body momentum distribution
//		functions based on intial state codes and multiplicity.
//
// NOTE:  Separate multiplicity-3 generators are not used, per V.Uzhinsky
//
// 20130619  Change singleton instance to be thread-local, to avoid collisions.
// 20130620  Address Coverity warnings about missing copy actions

#ifndef G4MultiBodyMomentumDist_h
#define G4MultiBodyMomentumDist_h 1

#include "globals.hh"

class G4VMultiBodyMomDst;
class G4NuclNucl3BodyMomDst;
class G4NuclNucl4BodyMomDst;
class G4HadNucl3BodyMomDst;
class G4HadNucl4BodyMomDst;


class G4MultiBodyMomentumDist {
public:
  ~G4MultiBodyMomentumDist();

  static const G4MultiBodyMomentumDist* GetInstance();

  // Return appropriate generator for initial state and multiplicity
  static const G4VMultiBodyMomDst* GetDist(G4int is, G4int mult) {
    return GetInstance()->ChooseDist(is, mult);
  }

  // Pass verbosity through to owned objects
  static void setVerboseLevel(G4int vb=0);

private:
  // Constructor is private for singleton
  G4MultiBodyMomentumDist();
  const G4VMultiBodyMomDst* ChooseDist(G4int is, G4int mult) const;

  void passVerbose(G4int verbose);	// Pass verbosity through instance

  static G4ThreadLocal G4MultiBodyMomentumDist* theInstance;	// Per thread

  // Generators for various initial/final state combinations
  G4NuclNucl3BodyMomDst* nn3BodyDst;	// N N to X Y Z
  G4NuclNucl4BodyMomDst* nn4BodyDst;    // N N to X Y Z W ...
  G4HadNucl3BodyMomDst*  hn3BodyDst;	// pi,K,Y,gamma N to X Y Z
  G4HadNucl4BodyMomDst*  hn4BodyDst;	// pi,K,Y,gamma N to X Y Z W ...

private:
  // Copying of modules is forbidden
  G4MultiBodyMomentumDist(const G4MultiBodyMomentumDist&);
  G4MultiBodyMomentumDist& operator=(const G4MultiBodyMomentumDist&);
};

#endif	/* G4MultiBodyMomentumDist_h */
