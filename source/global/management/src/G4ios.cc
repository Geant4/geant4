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
// $Id$
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

G4ThreadLocal G4strstreambuf *G4coutbuf_G4MT_TLS_ = 0;
G4ThreadLocal G4strstreambuf *G4cerrbuf_G4MT_TLS_ = 0;
G4ThreadLocal std::ostream *G4cout_G4MT_TLS_ = 0;
G4ThreadLocal std::ostream *G4cerr_G4MT_TLS_ = 0;
#define G4coutbuf (*G4coutbuf_G4MT_TLS_)
#define G4cerrbuf (*G4cerrbuf_G4MT_TLS_)
#define G4cout (*G4cout_G4MT_TLS_)
#define G4cerr (*G4cerr_G4MT_TLS_)

void G4iosInitialization()
{
  if (G4coutbuf_G4MT_TLS_ == 0) G4coutbuf_G4MT_TLS_ = new G4strstreambuf;
  if (G4cerrbuf_G4MT_TLS_ == 0) G4cerrbuf_G4MT_TLS_ = new G4strstreambuf;
  if (G4cout_G4MT_TLS_ == 0) G4cout_G4MT_TLS_ = new std::ostream(G4coutbuf_G4MT_TLS_);
  if (G4cerr_G4MT_TLS_ == 0) G4cerr_G4MT_TLS_ = new std::ostream(G4cerrbuf_G4MT_TLS_);
}

void G4iosFinalization()
{
  delete G4cout_G4MT_TLS_;
  delete G4cerr_G4MT_TLS_; 
  delete G4coutbuf_G4MT_TLS_;
  delete G4cerrbuf_G4MT_TLS_;
}

void __attribute__ ((constructor)) my_init(void)
{
  G4iosInitialization();
  //printf("G4ios is initialized\n");
}

void __attribute__ ((destructor)) my_fini(void)
{
  G4iosFinalization();
  //printf("G4ios is finalized\n");
}
