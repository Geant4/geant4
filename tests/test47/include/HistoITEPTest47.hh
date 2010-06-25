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
#ifndef HistoITEPTest47_H
#define HistoITEPTest47_H

#include "HistoTest47.hh"
#include "HistoEPTest47.hh"

#include "TFile.h"
#include "TH1F.h"
#include <string>

class G4VParticleChange;

class HistoITEPTest47 : public HistoTest47 {

public:

  HistoITEPTest47(std::string namePart, std::string nameMat, G4double momentum,
		  std::string nameGen);
  virtual ~HistoITEPTest47();

  virtual void fill(G4VParticleChange*, G4LorentzVector);
  virtual void write(G4double cross_sec, G4int nevt);

private:

  void initialize();
  void book();

  std::string           fileName;
  char                  tag1Name[60], tag2Name[24], tag3Name[40];
  double                dtheta, de;
  HistoEPTest47         epTest;
  std::vector<TH1F*>    hiKE11, hiKE12, hiCT11, hiCT12;
  std::vector<TH1F*>    hiKE21, hiKE22, hiCT21, hiCT22;
  std::vector<G4double> energies, emin, emax, angles, cthmin, cthmax, dcth;
  
};

#endif
