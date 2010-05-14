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

#include "G4CascadeChannel.hh"

std::vector<G4int> 
G4CascadeChannel::getQnums(G4int type)
{
  G4int bary=0, str=0, ch=0;
  std::vector<G4int> Qnums(3);

    switch(type) {
    case 1: // proton
      bary = 1;
      str = 0;
      ch = 1;
      break;
    case 2: // neutron
      bary = 1;
      str = 0;
      ch = 0;
      break;
    case 3: // pi+
      bary = 0;
      str = 0;
      ch = 1;
      break;
    case 5: // pi-
      bary = 0;
      str = 0;
      ch = -1;
      break;
    case 7: // pi0
      bary = 0;
      str = 0;
      ch = 0;
      break;
    case 11: // k+
      bary = 0;
      str = 1;
      ch = 1;
      break;
    case 13: // k-
      bary = 0;
      str = -1;
      ch = -1;
      break;
    case 15: // k0
      bary = 0;
      str = 1;
      ch = 0;
      break;
    case 17: // k0bar
      bary = 0;
      str = -1;
      ch = 0;
      break;
    case 21: // lambda
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 23: // sigma+
      bary = 1;
      str = -1;
      ch = 1;
      break;
    case 25: // sigma0
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 27: // sigma-
      bary = 1;
      str = -1;
      ch = -1;
      break;
    case 29: // xi0
      bary = 1;
      str = -2;
      ch = 0;
      break;
    case 31: // xi-
      bary = 1;
      str = -2;
      ch = -1;
      break;
    default:
      G4cout << " Unknown particle type " << type << G4endl;
    };

    Qnums[0] = bary;
    Qnums[1] = str;
    Qnums[2] = ch;
    return Qnums;
}
