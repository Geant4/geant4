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
// $Id: G4HtmlPPReporter.hh,v 1.1 2003/09/21 19:38:49 kurasige Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 
// ---------------------------------------------------------------
#ifndef G4HtmlPPReporter_h
#define G4HtmlPPReporter_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

#include "G4VParticlePropertyReporter.hh"

class G4HtmlPPReporter: public G4VParticlePropertyReporter
{
 public:
  //constructors
  G4HtmlPPReporter();
  
  //destructor
  virtual ~G4HtmlPPReporter();

 public:
  virtual void Print(const G4String& option="");
  
 private:
  void SparseOption(const G4String& option);
  void GenerateIndex();
  void GeneratePropertyTable(G4ParticleDefinition* );
  
  void PrintHeader(std::ofstream& );
  void PrintFooter(std::ofstream& );
  
  
 private:
  static const char *sTABLE, *eTABLE;
  static const char *sTR, *eTR;
  static const char *sTD, *eTD;
  static const char *sB, *eB;
  static const char *sLFONT, *eLFONT;
  static const char *sSYMBOL, *eSYMBOL;
  static const char *sSUP, *eSUP;
  static const char *sSUB, *eSUB;
  
 private:
  G4String  baseDir;
  G4String  comment;
};


#endif












