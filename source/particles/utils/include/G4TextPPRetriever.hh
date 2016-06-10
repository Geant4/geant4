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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TextPPRetriever.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ---------------------------------------------------------------
#ifndef G4TextPPRetriever_h
#define G4TextPPRetriever_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>
#include <fstream>       

#include "G4VParticlePropertyRetriever.hh"

class G4TextPPRetriever: public G4VParticlePropertyRetriever
{
 public:
  //constructors
  G4TextPPRetriever();
  
  //destructor
  virtual ~G4TextPPRetriever();
  
 public:
  virtual void Retrieve(const G4String& option="");

 protected:
  void SparseOption(const G4String& option);
  G4bool ModifyPropertyTable(const G4ParticleDefinition* );

 protected:
  G4String  baseDir;

};


#endif
