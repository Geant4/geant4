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
/// \file hadronic/Hadr02/include/UrQMDAntiBarionBuilder.hh
/// \brief Definition of the UrQMDAntiBarionBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   UrQMDAntiBarionBuilder
//
// Author: 2012  Andrea Dotti
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef UrQMDAntiBarionBuilder_h
#define UrQMDAntiBarionBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4VAntiBarionBuilder.hh"

#include "G4UrQMD1_3Model.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4VComponentCrossSection.hh"

class UrQMDAntiBarionBuilder : public G4VAntiBarionBuilder
{
public: 

  UrQMDAntiBarionBuilder();
  virtual ~UrQMDAntiBarionBuilder();

  virtual void Build(G4HadronElasticProcess * aP);
  virtual void Build(G4AntiProtonInelasticProcess * aP);
  virtual void Build(G4AntiNeutronInelasticProcess * aP);
  virtual void Build(G4AntiDeuteronInelasticProcess * aP);
  virtual void Build(G4AntiTritonInelasticProcess * aP);
  virtual void Build(G4AntiHe3InelasticProcess * aP);
  virtual void Build(G4AntiAlphaInelasticProcess * aP);
    
  inline void SetMinEnergy(G4double val) {fMin = val;}
  inline void SetMaxEnergy(G4double val) {fMax = val;}

private:

  G4UrQMD1_3Model * fModel; 
  G4VCrossSectionDataSet* fAntiNucleonData;
  G4VComponentCrossSection * fAntiNucleonXS;
  G4double fMin;
  G4double fMax;

};
#endif

