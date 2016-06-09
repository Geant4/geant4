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
// $Id: Analysis.hh,v 1.1 2007-11-16 14:29:32 kmura Exp $
// $Name: not supported by cvs2svn $
//
// ====================================================================
//   Analysis.hh
//
//                                         2007 Q
// ====================================================================
#ifndef ANALYSIS_MANAGER_H
#define ANALYSIS_MANAGER_H

#include "G4ThreeVector.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class TH1D;
class TH2D;

class Analysis {
private:
  static Analysis* myanalysis;

  TH2D* incident_map;
  TH1D* incident_x_hist;

  TH2D* dose_map;
  TH1D* dose_hist;

  Analysis(); // protected

  G4bool incidentFlag;

public:
  ~Analysis();

  static Analysis* GetAnalysis();

  void Update();
  void Clear();
  void Save(const G4String& fname);

  void FillIncident(const G4ThreeVector& p);
  void FillDose(const G4ThreeVector& p, G4double dedx);

  void ClearIncidentFlag();

};

// ====================================================================
//   inline functions
// ====================================================================

inline void Analysis::ClearIncidentFlag()
{ incidentFlag= false; }

#endif
