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
// G4coutDestination
//
// Class description:
//
// Cout/cerr buffer containers

// Authors: H.Yoshida, M.Nagamatu - November 1998
// --------------------------------------------------------------------
#ifndef G4COUTDESTINATION_HH
#define G4COUTDESTINATION_HH 1

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "globals.hh"

class G4coutDestination
{
 public:

  G4coutDestination() = default;
  virtual ~G4coutDestination() = default;

  // The type of the functions defining a transformation of the message.
  // The function manipulates the input message, for example, to add a
  // prefix:
  //    G4coutDestination::AddCoutTransformer(
  //       [](G4String& msg) -> G4bool { msg="PREFIX "+msg; return true; }
  //    );
  // Function should return false if message should not be processed
  // anymore and discarded
  //
  using Transformer = std::function<G4bool(G4String&)>;
  void AddCoutTransformer(const Transformer& t)
  {
    transformersCout.push_back(t);
  }
  void AddCoutTransformer(Transformer&& t) { transformersCout.push_back(t); }
  void AddCerrTransformer(const Transformer& t)
  {
    transformersCerr.push_back(t);
  }
  void AddCerrTransformer(Transformer&& t) { transformersCerr.push_back(t); }
  virtual void ResetTransformers();

  virtual G4int ReceiveG4cout(const G4String& msg);
  virtual G4int ReceiveG4cerr(const G4String& msg);
  // Derived class implements here handling of message.
  // For example, streaming on std::cout or file.
  // Return 0 for success, -1 otherwise

  G4int ReceiveG4cout_(const G4String& msg);
  // Method called by G4strbuf when need to handle a message

  G4int ReceiveG4cerr_(const G4String& msg);
  // Transformers cannot remove an error message from stream

 protected:

  G4MTGLOB_DLL static G4coutDestination* masterG4coutDestination;
  // For MT: if master G4coutDestination derived class wants to
  // intercept the thread outputs, derived class should set this pointer.
  // Needed for some G4UIsession like GUIs

  std::vector<Transformer> transformersCout;
  std::vector<Transformer> transformersCerr;
};

#endif
