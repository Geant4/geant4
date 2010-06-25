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
// $Id: G4CascadeChannel.hh,v 1.7 2010-06-25 09:41:52 gunter Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100514  M. Kelsey -- All functionality removed except quantum-number
//		validation functions.

#ifndef G4_CASCADE_CHANNEL_HH
#define G4_CASCADE_CHANNEL_HH

#include "globals.hh"
#include "G4FastVector.hh"
#include "G4ReactionProduct.hh"
#include <vector>

namespace G4CascadeChannel {
  std::vector<G4int> getQnums(G4int type);

  void CheckQnums(const G4FastVector<G4ReactionProduct,256> &vec,
		  G4int &vecLen,
		  G4ReactionProduct &currentParticle,
		  G4ReactionProduct &targetParticle,
		  G4double Q, G4double B, G4double S);
}

#endif
