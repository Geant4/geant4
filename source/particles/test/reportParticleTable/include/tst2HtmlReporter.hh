// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2HtmlReporter.hh,v 1.1 1999-06-17 04:41:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2HtmlReporter_h
#define tst2HtmlReporter_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>

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

    void PrintHeader(ofstream& );
    void PrintFooter(ofstream& );


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
