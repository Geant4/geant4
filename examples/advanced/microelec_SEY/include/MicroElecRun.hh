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
/// \file MicroElecRun.cc
/// \brief Implementation of the MicroElecRun class
//
//

#ifndef MicroElecRun_h
#define MicroElecRun_h 1

#include "globals.hh"
#include "G4Run.hh"

#include "G4THitsMap.hh"

class G4Event;

class MicroElecRun : public G4Run
{
public:

  MicroElecRun();
  ~MicroElecRun() override;
  
  void RecordEvent(const G4Event*) override;
  void Merge(const G4Run*) override;
  
  G4double GetElecPrimScorer()     { return ElecPrimScorer; }
  G4double GetElecSecoScorer()     { return ElecSecoScorer; }
  G4double GetElecSup50Scorer() { return ElecSup50Scorer; }
  G4double GetElecTotaScorer()     { return ElecTotaScorer; }
  G4double GetElecEneIncPart() { return ElecEneIncPart; }
  
  void SetElecPrimScorer(G4double elecprimscorer) { ElecPrimScorer= elecprimscorer; }
  void SetElecSecoScorer(G4double elecsecoscorer) { ElecSecoScorer= elecsecoscorer; }
  void SetElecSup50Scorer(G4double elecsup50scorer) { ElecSup50Scorer= elecsup50scorer; }
  void SetElecTotaScorer(G4double electotascorer) { ElecTotaScorer= electotascorer; }
  void SetElecEneIncPart(G4double eleceneincPart) { ElecEneIncPart = eleceneincPart; }

  void AddElecPrimScorer(G4double elecprimscorer) { ElecPrimScorer += elecprimscorer; }
  void AddElecSecoScorer(G4double elecsecoscorer) { ElecSecoScorer += elecsecoscorer; }
  void AddElecSup50Scorer(G4double elecsup50scorer) { ElecSup50Scorer += elecsup50scorer; }
  void AddElecTotaScorer(G4double electotascorer) { ElecTotaScorer += electotascorer; }



private:
  G4double ElecPrimScorer;
  G4double ElecSecoScorer;
  G4double ElecTotaScorer;
  G4double ElecSup50Scorer;
  G4double ElecEneIncPart;
	
  
};

#endif

