//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef Tst69INCLXXTallyAnalysis_hh
#define Tst69INCLXXTallyAnalysis_hh

#include "G4INCLXXVInterfaceTally.hh"
#include "G4String.hh"

class Tst69INCLXXTallyAnalysis : public G4INCLXXVInterfaceTally {
  public:
    Tst69INCLXXTallyAnalysis(G4String const &physList, const G4bool histos=true, const G4bool ntuple=false);
    virtual ~Tst69INCLXXTallyAnalysis();
    virtual void Open();
    virtual void Close();
    virtual void Tally(G4HadProjectile const &aTrack, G4Nucleus const &theNucleus, G4HadFinalState const &result);

  protected:
    G4int whichIncomingParticle(const G4int A, const G4int Z);
    G4int whichOutgoingParticle(const G4int A, const G4int Z);
    G4int whichEnergySlice(const G4double E);

    G4String filenameStem;
    G4String filename;

    G4bool factoryOn;
    G4bool withHistos;
    G4bool withNtuple;

    G4int ntupleID;
    G4int colIDAp;
    G4int colIDZp;
    G4int colIDEp;
    G4int colIDAt;
    G4int colIDZt;

    static const int nIncPart = 4;
    static const int nOutPart = 5;
    static const int nESlice = 3;
    static const char *incPartDict[nIncPart];
    static const char *outPartDict[nOutPart];
    static const char *incPartLongDict[nIncPart];
    static const char *outPartLongDict[nOutPart];
    static const char *eSliceDict[nESlice];
    static const float eSliceLow[nESlice+1];
    static const int nBins;
    static const int nZBins;

    G4int hID_ioe[nIncPart][nOutPart][nESlice];
    G4int hID_Z_ie[nIncPart][nESlice];
    G4int hID_i[nIncPart];
};

#endif
