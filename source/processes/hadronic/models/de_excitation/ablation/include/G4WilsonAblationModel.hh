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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4WilsonAblationModel_h
#define G4WilsonAblationModel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4WilsonAblationModel.hh
//
// Version:		1.0
// Date:		08/12/2009
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// 08 December 2009, P R Truscott, QinetiQ Ltd, UK
// Ver 1.0
// Introduced vector of evaporation channels and evaporation factory.  Also
// copied directly over the SumProbabilities class in G4Evaporation.hh at
// version 9.2.r9, just in cases there's any subtle differences.  See .cc
// file comments to see impact of the rest of the changes.
//
// 04 October 2014, D Mancusi
// Moved theChannels and theChannelFactory to the base class, since they seem
// to be common to all classes derived from G4VEvaporation.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VEvaporation.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4EvaporationFactory.hh"


////////////////////////////////////////////////////////////////////////////////
//
class G4WilsonAblationModel : public G4VEvaporation
{
  public:
    G4WilsonAblationModel();
    virtual ~G4WilsonAblationModel();
    
    typedef std::vector<G4ParticleDefinition*> VectorOfFragmentTypes;

    G4FragmentVector * BreakItUp (const G4Fragment &theNucleus);
    void SetProduceSecondaries (G4bool);
    G4bool GetProduceSecondaries ();
    void SetVerboseLevel (G4int);
    G4int GetVerboseLevel ();

  private:
    void SelectSecondariesByEvaporation (G4Fragment*);
    void SelectSecondariesByDefault (G4ThreeVector);
    void PrintWelcomeMessage ();

  private:
    G4bool                 produceSecondaries;
    G4int                  verboseLevel;
    G4double               B;
    G4int                  nFragTypes;
    G4ParticleDefinition  *fragType[6];
    G4FragmentVector      *fragmentVector;
    VectorOfFragmentTypes  evapType;
    G4double               fSig[200];
};
////////////////////////////////////////////////////////////////////////////////
//
inline void G4WilsonAblationModel::SetProduceSecondaries 
  (G4bool produceSecondaries1)
  {produceSecondaries = produceSecondaries1;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4bool G4WilsonAblationModel::GetProduceSecondaries ()
  {return produceSecondaries;}
////////////////////////////////////////////////////////////////////////////////
//
inline void G4WilsonAblationModel::SetVerboseLevel (G4int verboseLevel1)
  {verboseLevel = verboseLevel1;}
////////////////////////////////////////////////////////////////////////////////
//
inline G4int G4WilsonAblationModel::GetVerboseLevel ()
  {return verboseLevel;}
////////////////////////////////////////////////////////////////////////////////
//
#endif
