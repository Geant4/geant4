#ifndef G4WilsonAblationModel_h
#define G4WilsonAblationModel_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4WilsonAblationModel.hh
//
// Version:		B.1
// Date:		15/04/04
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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation  is  the  intellectual   property  of *
// * QinetiQ.   Rights  to use  this  software  are  granted  to  the *
// * European Space Agency  under  the  provisions  of  Clause 39  of *
// * the    Standard    Conditions   for  Contract.   ESA    contract *
// * 17191/03/NL/LvH  provides  the  rights to distribute source code *
// * through  the  Geant4 Collaboration Agreement for general  public *
// * use.  Some elements of this source code may be derived from code *
// * developed   by   other   member   institutes   of   the   Geant4 *
// * Collaboration, and the provision of  their code under the Geant4 *
// * Collaboration is acknowledged.                                   *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include <vector>
////////////////////////////////////////////////////////////////////////////////
//
class G4WilsonAblationModel : public G4VEvaporation
{
  public:
    G4WilsonAblationModel();
    ~G4WilsonAblationModel();
    
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

    class SumProbabilities : 
      public std::binary_function<G4double,G4double,G4double>
    {
      public:
        SumProbabilities() : total(0.0) {}
        G4double operator() (G4double& /* probSoFar */, G4VEvaporationChannel*& frag)
        {
          total += frag->GetEmissionProbability();
          return total;
        }

        G4double GetTotal() { return total; }
      public:
      G4double total;
  };
                                                                                


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
