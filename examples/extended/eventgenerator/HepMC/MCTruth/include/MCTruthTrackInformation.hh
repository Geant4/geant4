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
/// \file eventgenerator/HepMC/MCTruth/include/MCTruthTrackInformation.hh
/// \brief Definition of the MCTruthTrackInformation class
//
//
// $Id: MCTruthTrackInformation.hh 103182 2017-03-21 10:36:09Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 - MCTruthTrackInformation class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
// Date  : 2006-03-06
//
// --------------------------------------------------------------
#ifndef INCLUDE_MCTRUTHTRACKINFORMATION_H 
#define INCLUDE_MCTRUTHTRACKINFORMATION_H 1

#include "G4Types.hh"
#include "G4VUserTrackInformation.hh"

class MCTruthTrackInformation : public G4VUserTrackInformation
{
public:

  MCTruthTrackInformation(); 

  virtual ~MCTruthTrackInformation();

  void Print() const {};

  void SetDirectParent(G4bool value) { fDirectParent = value; }
  G4bool GetDirectParent() const { return fDirectParent; }

private:

  G4bool fDirectParent;

};

#endif // INCLUDE_MCTRUTHTRACKINFORMATION_H
