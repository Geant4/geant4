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
// $Id: tst2HtmlReporter.hh,v 1.3 2001-07-11 10:02:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2HtmlReporter_h
#define tst2HtmlReporter_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"

#include "tst2VParticleReporter.hh"

class tst2HtmlReporter: public tst2VParticleReporter
{
 public:
  //constructors
    tst2HtmlReporter();

  //destructor
    virtual ~tst2HtmlReporter();

 public:
	virtual void Print(const tst2ParticleContainer& container, 
                       const G4String& option="");

 private:
    void SparseOption(const G4String& option);
    void GenerateIndex();
    void GeneratePropertyTable(G4ParticleDefinition* );

    void PrintHeader(G4std::ofstream& );
    void PrintFooter(G4std::ofstream& );


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
