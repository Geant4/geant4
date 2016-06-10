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
// $Id: G4InuclParticleNames.hh 69638 2013-05-09 04:26:00Z mkelsey $
//
// Defines enums to map G4InuclElementaryParticle type codes to human
// readable names.  Meant to replace similar local enums scattered through
// the code.

#include "G4InuclParticleNames.hh"

using namespace G4InuclParticleNames;


// Convert enum entries to strings

const char* G4InuclParticleNames::nameLong(G4int ptype) {
  switch (ptype) {
  case nuclei: return "nuclei"; break;
  case proton: return "proton"; break;
  case neutron: return "neutron"; break;
  case pionPlus: return "pionPlus"; break;
  case pionMinus: return "pionMinus"; break;
  case pionZero: return "pionZero"; break;
  case photon: return "photon"; break;
  case kaonPlus: return "kaonPlus"; break;
  case kaonMinus: return "kaonMinus"; break;
  case kaonZero: return "kaonZero"; break;
  case kaonZeroBar: return "kaonZeroBar"; break;
  case lambda: return "lambda"; break;
  case sigmaPlus: return "sigmaPlus"; break;
  case sigmaZero: return "sigmaZero"; break;
  case sigmaMinus: return "sigmaMinus"; break;
  case xiZero: return "xiZero"; break;
  case xiMinus: return "xiMinus"; break;
  case omegaMinus: return "omegaMinus"; break;
  case deuteron: return "deuteron"; break;
  case triton: return "triton"; break;
  case He3: return "He3"; break;
  case alpha: return "alpha"; break;
  case antiProton: return "antiProton"; break;
  case antiNeutron: return "antiNeutron"; break;
  case antiDeuteron: return "antiDeuteron"; break;
  case antiTriton: return "antiTriton"; break;
  case antiHe3: return "antiHe3"; break;
  case antiAlpha: return "antiAlpha"; break;
  case diproton: return "diproton"; break;
  case unboundPN: return "unboundPN"; break;
  case dineutron: return "dineutron"; break;
  case electronNu: return "electronNu"; break;
  case muonNu: return "muonNu"; break;
  case tauNu: return "tauNu"; break;
  case antiElectronNu: return "antiElectronNu"; break;
  case antiMuonNu: return "antiMuonNu"; break;
  case antiTauNu: return "antiTauNu"; break;
  case WMinus: return "WMinus"; break;
  case WPlus: return "WPlus"; break;
  case Zzero: return "Zzero"; break;
  case electron: return "electron"; break;
  case muonMinus: return "muonMinus"; break;
  case tauMinus: return "tauMinus"; break;
  case positron: return "positron"; break;
  case muonPlus: return "muonPlus"; break;
  case tauPlus: return "tauPlus"; break;
  default: ;
  }
  return "UNKNOWN";
}

const char* G4InuclParticleNames::nameShort(G4int ptype) {
  switch (ptype) {
  case nuc: return "nuc"; break;
  case pro: return "pro"; break;
  case neu: return "neu"; break;
  case pip: return "pip"; break;
  case pim: return "pim"; break;
  case pi0: return "pi0"; break;
  case gam: return "gam"; break;
  case kpl: return "kpl"; break;
  case kmi: return "kmi"; break;
  case k0: return "k0"; break;
  case k0b: return "k0b"; break;
  case lam: return "lam"; break;
  case sp: return "sp"; break;
  case s0: return "s0"; break;
  case sm: return "sm"; break;
  case xi0: return "xi0"; break;
  case xim: return "xim"; break;
  case om: return "om"; break;
  case deu: return "deu"; break;
  case ap: return "ap"; break;
  case an: return "an"; break;
  case ade: return "ade"; break;
  case atr: return "atr"; break;
  case ahe: return "ahe"; break;
  case aal: return "aal"; break;
  case pp: return "pp"; break;
  case pn: return "pn"; break;
  case nn: return "nn"; break;
  case enu: return "enu"; break;
  case mnu: return "mnu"; break;
  case tnu: return "tnu"; break;
  case aenu: return "aenu"; break;
  case amnu: return "amnu"; break;
  case atnu: return "atnu"; break;
  case wm: return "wm"; break;
  case wp: return "wp"; break;
  case z0: return "z0"; break;
  case ele: return "ele"; break;
  case mum: return "mum"; break;
  case tm: return "tm"; break;
  case pos: return "pos"; break;
  case mup: return "mup"; break;
  case tp: return "tp"; break;
  default: ;
  }
  return "?";
}
