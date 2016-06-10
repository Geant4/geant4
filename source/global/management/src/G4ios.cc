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
// $Id: G4ios.cc 70021 2013-05-22 07:55:29Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
// G4ios.cc
//
// History 1998 Nov. 3 Masayasu Nagamatu

#include "G4ios.hh"
#include "G4strstreambuf.hh"

#ifdef G4MULTITHREADED
  G4ThreadLocal G4strstreambuf *G4coutbuf_p = 0;
  G4ThreadLocal G4strstreambuf *G4cerrbuf_p = 0;
  G4ThreadLocal std::ostream *G4cout_p = 0;
  G4ThreadLocal std::ostream *G4cerr_p = 0;
  #define G4coutbuf (*G4coutbuf_p)
  #define G4cerrbuf (*G4cerrbuf_p)
  #define G4cout (*G4cout_p)
  #define G4cerr (*G4cerr_p)

  void G4iosInitialization()
  {
    if (G4coutbuf_p == 0) G4coutbuf_p = new G4strstreambuf;
    if (G4cerrbuf_p == 0) G4cerrbuf_p = new G4strstreambuf;
    if (G4cout_p == 0) G4cout_p = new std::ostream(G4coutbuf_p);
    if (G4cerr_p == 0) G4cerr_p = new std::ostream(G4cerrbuf_p);
  }

  void G4iosFinalization()
  {
    delete G4cout_p;
    delete G4cerr_p; 
    delete G4coutbuf_p;
    delete G4cerrbuf_p;
  }

#else  // Sequential

  G4strstreambuf G4coutbuf;
  G4strstreambuf G4cerrbuf;
  std::ostream G4cout(&G4coutbuf);
  std::ostream G4cerr(&G4cerrbuf);

  void G4iosInitialization() {}
  void G4iosFinalization() {}

#endif
