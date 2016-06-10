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
// $Id: G4TypeKeyT.hh 69802 2013-05-15 14:52:57Z gcosmo $
//
// Type key class
//
// Jane Tinslay, September 2006
//
#ifndef G4TYPEKEYT_HH
#define G4TYPEKEYT_HH

#include "G4TypeKey.hh"

template <typename T>
class G4TypeKeyT : public G4TypeKey { 

public:
  
  G4TypeKeyT() {
    static G4ThreadLocal Key *pkey = 0 ;
    if (!pkey) { pkey = new Key; *pkey = NextKey(); }
    Key &key = *pkey;
    fMyKey = key;
  }

  virtual ~G4TypeKeyT() {}

};

#endif
