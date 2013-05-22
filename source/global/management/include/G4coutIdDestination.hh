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
// ********************************************************************//
//
// $Id: G4coutDestination.hh 66241 2014-05-16 00:00:42Z adotti $
//
//
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Prepends to G4cout and G4cerr an id
// Used by default for Geant4MT
// Optionally it buffers all output and sends it to
// output on request.
// By default streaming is on std::cout and std::cerr
// but this can be changed.
// Warning: this option stores streams in memory.
//
// Example:
//  #include "G4coutIdDestination.hh"
//  #include "G4ios.hh"
//  int main(int,char**) {
//      G4UImanager::GetUIpointer();
//      G4iosInitialization();
//      G4int id = 1;
//      G4coutIdDestination myout(id);
#ifndef G4COUTIDDESTINATION_HH
#define G4COUTIDDESTINATION_HH

#include "G4coutDestination.hh"
#include <iostream>
#include <sstream>

class G4coutIdDestination : public G4coutDestination
{
public:
    G4coutIdDestination( const G4int& id, std::ostream& cout=std::cout, std::ostream&  cerr=std::cerr );
    virtual ~G4coutIdDestination();
    virtual G4int ReceiveG4cout(const G4String&);
    virtual G4int ReceiveG4cerr(const G4String&);
    void EnableBuffering(G4bool flag=true);
    void DumpBuffer();
private:
    std::ostream& finalcout;
    std::ostream& finalcerr;
    const G4int id;
    G4bool useBuffer;
    std::ostringstream cout_buffer;
    std::ostringstream cerr_buffer;
};

#endif // G4COUTIDDESTINATION_HH
