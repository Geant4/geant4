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
// $Id: G4PiKBuilder.hh,v 1.2 2005/11/25 15:38:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
//
//----------------------------------------------------------------------------
//
#ifndef G4PiKBuilder_h
#define G4PiKBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include <vector>

class G4PiKBuilder
{
  public: 
    G4PiKBuilder();
    virtual ~G4PiKBuilder();

  public: 
    void Build();
    void RegisterMe(G4VPiKBuilder * aB) {theModelCollections.push_back(aB);}

  private:
    G4HadronElasticProcess* thePionPlusElasticProcess;
    G4HadronElasticProcess* thePionMinusElasticProcess;
    G4HadronElasticProcess* theKaonPlusElasticProcess;
    G4HadronElasticProcess* theKaonMinusElasticProcess;
    G4HadronElasticProcess* theKaonZeroLElasticProcess;
    G4HadronElasticProcess* theKaonZeroSElasticProcess;

    G4PionPlusInelasticProcess*  thePionPlusInelastic;
    G4PionMinusInelasticProcess* thePionMinusInelastic;
    G4KaonPlusInelasticProcess*  theKaonPlusInelastic;
    G4KaonMinusInelasticProcess* theKaonMinusInelastic;
    G4KaonZeroLInelasticProcess* theKaonZeroLInelastic;
    G4KaonZeroSInelasticProcess* theKaonZeroSInelastic;
     
    std::vector<G4VPiKBuilder *> theModelCollections;

    G4bool wasActivated;
};

// 2002 by J.P. Wellisch

#endif

