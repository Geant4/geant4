// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gspart.cc,v 1.2 1999-11-15 10:39:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ProcessManager.hh"
#include "G3toG4.hh"
#include "G3PartTable.hh"

void PG4gspart(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgspart);

    // interpret the parameters
    G4String chnpar = Spar[0];
    G4int ipart = Ipar[0];
    G4int itrtyp = Ipar[1];
    G4int nwb = Ipar[2];
    G4double amass = Rpar[0];
    G4double charge = Rpar[1];
    G4double tlife = Rpar[2];
    G4double *ubuf = &Rpar[3];

    G4gspart(ipart,chnpar,itrtyp,amass,charge,tlife,ubuf,nwb);
}

void G4gspart(G4int ipart, G4String chnpar, G4int itrtyp, G4double amass,
              G4double charge, G4double tlife, G4double*, G4int)
{
}

#ifdef OMIT_CODE
void G4gspart(G4int ipart, G4String chnpar, G4int itrtyp, G4double amass,
              G4double charge, G4double tlife, G4double*, G4int)
{
    // Handle conversion of itrtyp into an appropriate ProcessManager
    G4ProcessManager *mgr = NULL;
    switch (itrtyp) {
    case 1:
        // gamma
// $$$        mgr = gammaProcessManager;
        break;
    case 2:
        // electron
// $$$        mgr = electronProcessManager;
        break;
    case 3:
        // neutron
// $$$        mgr = neutronProcessManager;
        break;
    case 4:
        // hadron
// $$$        mgr = hadronProcessManager;
        break;
    case 5:
        // muon
// $$$        mgr = muonProcessManager;
        break;
    case 6:
        // geantino
// $$$        mgr = geantinoProcessManager;
        break;
    case 7:
        // heavy ion
// $$$        mgr = heavyIonProcessManager;
        break;
    case 8:
        // light ion
// $$$        mgr = lightIonProcessManager;
        break;
    default:
// $$$        mgr = theProcessManager;
        break;
    }

    // Create the particle; provide parameters and a process mgr.
    G4ParticleDefinition *part = new G4ParticleDefinition(chnpar);
    part->SetProcessManager(mgr);
    part->SetPDGMass(amass);
    part->SetPDGCharge(charge);
    part->SetPDGLifeTime(tlife);

    // add to particle table
    G3Part.put(&ipart, part);
}
#endif
