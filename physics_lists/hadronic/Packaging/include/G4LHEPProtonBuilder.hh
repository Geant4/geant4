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
#ifndef G4LHEPProtonBuilder_h
#define G4LHEPProtonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4LElastic.hh"   
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"

class G4LHEPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4LHEPProtonBuilder();
    virtual ~G4LHEPProtonBuilder();

  public: 
    virtual void Build(G4ProtonInelasticProcess & aP);
    virtual void Build(G4HadronElasticProcess & aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin = aM;
    }

  private:
    G4LElastic * theElasticModel;
    G4LEProtonInelastic * theLEProtonModel;
    G4HEProtonInelastic * theHEProtonModel;
    
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

