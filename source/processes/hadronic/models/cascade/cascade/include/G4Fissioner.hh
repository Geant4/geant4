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
// $Id: G4Fissioner.hh 71954 2013-06-29 04:40:40Z mkelsey $
//
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100728  M. Kelsey -- Move G4FissionStore to data member and reuse
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20110801  M. Kelsey -- Pass C arrays to ::potentialMinimization()
// 20130129  M. Kelsey -- Put buffer for output nuclei here for thread-safety
// 20130620  Address Coverity complaint about missing copy actions
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment
// 20130628  Drop local list of fragments; add directly to output

#ifndef G4FISSIONER_HH
#define G4FISSIONER_HH

#include "G4CascadeDeexciteBase.hh"
#include "G4FissionStore.hh"
#include "G4InuclNuclei.hh"
#include <vector>


class G4Fissioner : public G4CascadeDeexciteBase {
public:
  G4Fissioner() : G4CascadeDeexciteBase("G4Fissioner") {;}
  virtual ~G4Fissioner() {;}

  virtual void deExcite(const G4Fragment& target, G4CollisionOutput& output);

private: 
  G4FissionStore fissionStore;

  G4double getC2(G4int A1, 
		 G4int A2, 
		 G4double X3, 
		 G4double X4, 
		 G4double R12) const; 

  G4double getZopt(G4int A1, 
		   G4int A2, 
		   G4int ZT, 
                   G4double X3, 
		   G4double X4, 
		   G4double R12) const;
		    
  void potentialMinimization(G4double& VP, 
			     G4double (&ED)[2], 
			     G4double& VC,
			     G4int AF, 
			     G4int AS, 
			     G4int ZF, 
			     G4int ZS,
			     G4double AL1[2], 
			     G4double BET1[2], 
			     G4double& R12) const; 

private:
  // Copying of modules is forbidden
  G4Fissioner(const G4Fissioner&);
  G4Fissioner& operator=(const G4Fissioner&);
};        

#endif /* G4FISSIONER_HH */
