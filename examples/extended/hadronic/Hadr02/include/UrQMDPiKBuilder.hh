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
/// \file hadronic/Hadr02/include/UrQMDPiKBuilder.hh
/// \brief Definition of the UrQMDPiKBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   UrQMDPiKBuilder
//
// Author: 2012 Andrea Dotti
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef UrQMDPiKBuilder_h
#define UrQMDPiKBuilder_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include "G4UrQMD1_3Model.hh"

class UrQMDPiKBuilder : public G4VPiKBuilder
{
public: 

  UrQMDPiKBuilder();
  virtual ~UrQMDPiKBuilder();
 
  virtual void Build(G4HadronElasticProcess * aP);
  virtual void Build(G4PionPlusInelasticProcess * aP);
  virtual void Build(G4PionMinusInelasticProcess * aP);
  virtual void Build(G4KaonPlusInelasticProcess * aP);
  virtual void Build(G4KaonMinusInelasticProcess * aP);
  virtual void Build(G4KaonZeroLInelasticProcess * aP);
  virtual void Build(G4KaonZeroSInelasticProcess * aP);
    
  inline void SetMinEnergy(G4double aM) {fMin = aM;}
  inline void SetMaxEnergy(G4double aM) {fMax = aM;}

private:

  G4UrQMD1_3Model * fModel;
  G4double fMin;
  G4double fMax;

};

#endif

