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
#ifndef G4PiKBuilder_h
#define G4PiKBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include "g4std/vector"

class G4PiKBuilder
{
  public: 
    G4PiKBuilder();
    virtual ~G4PiKBuilder();

  public: 
    void Build();
    void RegisterMe(G4VPiKBuilder * aB) {theModelCollections.push_back(aB);}

  private:
    G4HadronElasticProcess thePionPlusElasticProcess;
    G4HadronElasticProcess thePionMinusElasticProcess;
    G4HadronElasticProcess theKaonPlusElasticProcess;
    G4HadronElasticProcess theKaonMinusElasticProcess;
    G4HadronElasticProcess theKaonZeroLElasticProcess;
    G4HadronElasticProcess theKaonZeroSElasticProcess;

    G4PionPlusInelasticProcess  thePionPlusInelastic;
    G4PionMinusInelasticProcess thePionMinusInelastic;
    G4KaonPlusInelasticProcess  theKaonPlusInelastic;
    G4KaonMinusInelasticProcess theKaonMinusInelastic;
    G4KaonZeroLInelasticProcess theKaonZeroLInelastic;
    G4KaonZeroSInelasticProcess theKaonZeroSInelastic;
     
    G4std::vector<G4VPiKBuilder *> theModelCollections;

};

// 2002 by J.P. Wellisch

#endif

