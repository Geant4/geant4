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
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TextPPRetriever.hh,v 1.1 2004/03/11 09:47:44 kurasige Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
  G4bool ModifyPropertyTable(G4ParticleDefinition* );

 protected:
  G4String  baseDir;

};


#endif
