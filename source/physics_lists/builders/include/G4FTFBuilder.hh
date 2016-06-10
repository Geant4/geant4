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
// $Id: G4FTFBuilder.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4FTFBuilder
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4FTFBuilder_h
#define G4FTFBuilder_h 1

#include "globals.hh"
#include "G4VHadronModelBuilder.hh"

class G4FTFModel;
class G4ExcitedStringDecay;
class G4PreCompoundModel;
class G4GeneratorPrecompoundInterface;
class G4LundStringFragmentation;

class G4FTFBuilder : public G4VHadronModelBuilder
{
public: 

  G4FTFBuilder(const G4String& name ="FTF",
	       G4PreCompoundModel* p = 0); 

  virtual ~G4FTFBuilder();

protected:

  virtual G4HadronicInteraction* BuildModel();

private:

  // copy constructor and hide assignment operator
  G4FTFBuilder(G4FTFBuilder &);
  G4FTFBuilder & operator=(const G4FTFBuilder &right);

  G4FTFModel*            fStringModel;
  G4ExcitedStringDecay*  fStringDecay;
  G4PreCompoundModel*    fPreCompound;

  G4GeneratorPrecompoundInterface* fPrecoInterface;
  G4LundStringFragmentation*       fLund;
};

#endif

