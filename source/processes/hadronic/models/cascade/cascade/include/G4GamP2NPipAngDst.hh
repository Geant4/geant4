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
// $Id: G4GamP2NPipAngDst.hh 67433 2013-02-20 21:27:31Z mkelsey $
// Author:  Dennis Wright (SLAC)
// Date:    28 January 2013
//
// Description: class containing numerically integrated angular distributions
//              in the CM for the gamma p -> n pi+ reaction
//
// 20130219	M. Kelsey: Inherit from templated base, remove ctor arguments


#ifndef G4GamP2NPipAngDst_h
#define G4GamP2NPipAngDst_h 1

#include "G4NumIntTwoBodyAngDst.hh"


class G4GamP2NPipAngDst : public G4NumIntTwoBodyAngDst<15,19> {
public:
  G4GamP2NPipAngDst(G4int verbose=0);
  virtual ~G4GamP2NPipAngDst() {}
};

#endif
