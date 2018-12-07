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
//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4AtomicShells_XDB_EADL.hh                                        //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   9 August 2018                                                     //
//                                                                            //
//  Description: Class containing number of shells, electron configurations   //
//               and binding energies for atoms from Z = 1 to Z = 120.        //
//               Most entries are taken from the X-ray Data Book, with        //
//               unmeasured values supplied by Kibedi.  See documentation     //
//               file in G4EMLOW7.3 or later, directory fluor/                //  
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef G4ATOMICSHELLS_XDB_EADL_H
#define G4ATOMICSHELLS_XDB_EADL_H

#include "globals.hh"


class G4AtomicShells_XDB_EADL
{
  public :  //with description

    static G4int    GetNumberOfShells    (G4int Z);
    static G4int    GetNumberOfElectrons (G4int Z, G4int SubshellNb);
    static G4double GetBindingEnergy     (G4int Z, G4int SubshellNb);
    static G4double GetTotalBindingEnergy(G4int Z);

  private :

#ifdef G4VERBOSE
    static G4int PrintErrorZ(G4int Z, const G4String&);
    static G4int PrintErrorShell(G4int Z, G4int SubshellNb, const G4String&);
#endif

    G4AtomicShells_XDB_EADL(const G4AtomicShells_XDB_EADL&) = delete;
    const G4AtomicShells_XDB_EADL& operator=(const G4AtomicShells_XDB_EADL&) = delete;

    static const G4int    fNumberOfShells   [121];
    static const G4int    fIndexOfShells    [121];
    static const G4int    fNumberOfElectrons[2171];
    static const G4double fBindingEnergies  [2171];
};

#endif   // end of G4AtomicShells_XDB_EADL.hh

