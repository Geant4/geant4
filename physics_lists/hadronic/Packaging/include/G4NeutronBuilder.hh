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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#ifndef G4NeutronBuilder_h
#define G4NeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include <vector>

class G4NeutronBuilder
{
  public: 
    G4NeutronBuilder();
    virtual ~G4NeutronBuilder();

  public: 
    void Build();
    void RegisterMe(G4VNeutronBuilder * aB) {theModelCollections.push_back(aB);}

  private:
    G4NeutronInelasticProcess * theNeutronInelastic;
    G4HadronFissionProcess * theNeutronFission;
    G4HadronCaptureProcess  * theNeutronCapture;
    
    std::vector<G4VNeutronBuilder *> theModelCollections;

    G4bool wasActivated;

};

// 2002 by J.P. Wellisch

#endif

