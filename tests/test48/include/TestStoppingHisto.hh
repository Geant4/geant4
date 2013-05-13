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
#ifndef TestStoppingHisto_h
#define TestStoppingHisto_h

#include <string>
#include <vector>

#include "G4TrackVector.hh"

#include "TFile.h"
#include "TH1F.h"

// fwd declaration
class G4VParticleChange;
class G4DynamicParticle;

class TestStoppingHisto
{

public:
      
  // ctor & dtor
  //
  TestStoppingHisto( std::string beam, std::string target, std::string model) 
    : fJobID(-1), fBeam(beam), fTarget(target), fModel(model) { Init(); }
  ~TestStoppingHisto();
      
  void FillEvt( G4VParticleChange* );
  void FillEvtKaonBeam( G4VParticleChange* );
  void FillEvtMuonMinusBeam( G4TrackVector * );
  void FillEvtSigmaBeam( G4VParticleChange* );
  int  Topology( G4VParticleChange* );
  void Write( int ) ;
  void SetJobID( int id ) { fJobID=id; return; }

protected:
   int fJobID; 

private:
      
  void Init();
  void InitHistoGeneral();
  void InitPionMinus();
  void InitKaonMinus();
  void InitMuonMinus();
  void InitSigmaMinus();
      
  // void FillEvtAntiProton( const G4DynamicParticle* );
      
  // data members
  //
  std::string        fBeam;
  std::string        fTarget;
  std::string        fModel;
  std::string        fHistoTitle;
  std::string        ptag;
  std::vector<TH1F*> fHisto; 
  std::vector<TH1F*> fMuHisto; 

};

#endif
