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
// $Id: G4MuonMinusCapture.hh,v 1.23 2008-10-02 20:57:52 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4MuonMinusCapture
//
// Author V.Ivanchenko 25 April 2012 on base of G4MuonMinusCaptureAtRest
//
//
// Class Description:
//
// Stopping of mu-
//
// Modifications: 
//
//------------------------------------------------------------------------

#ifndef G4MuonMinusCapture_h
#define G4MuonMinusCapture_h 1
 
#include "globals.hh"
#include "G4HadronStoppingProcess.hh"
#include "G4ParticleDefinition.hh"

class G4VPreCompoundModel;

class G4MuonMinusCapture : public G4HadronStoppingProcess 
{ 
public:
 
  G4MuonMinusCapture(G4VPreCompoundModel* ptr=0);

  ~G4MuonMinusCapture();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void ProcessDescription(std::ostream& outFile) const;

private:

  // hide assignment operator as private 
  G4MuonMinusCapture& operator=(const G4MuonMinusCapture &right);
  G4MuonMinusCapture(const G4MuonMinusCapture& );

};

#endif
 





