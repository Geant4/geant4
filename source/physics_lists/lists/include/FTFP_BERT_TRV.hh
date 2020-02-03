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
//
//---------------------------------------------------------------------------
//
// ClassName:  FTFP_BERT_TRV
//
// Author: 2009  John Apostolakis
//   created from FTFP_BERT - changing threshold between BERT and FTF
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef FTFP_BERT_TRV_h
#define FTFP_BERT_TRV_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"


class FTFP_BERT_TRV: public G4VModularPhysicsList
{
public:
  FTFP_BERT_TRV(G4int ver = 1);
  virtual ~FTFP_BERT_TRV()=default;
  
  FTFP_BERT_TRV(const FTFP_BERT_TRV &) = delete;
  FTFP_BERT_TRV & operator=(const FTFP_BERT_TRV &)=delete;
};

#endif



