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
#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4ProtonInelasticCrossSection.hh"

class G4QGSCProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4QGSCProtonBuilder();
    virtual ~G4QGSCProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4QGSParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

