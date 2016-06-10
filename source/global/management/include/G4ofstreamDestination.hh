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
// $Id: G4coutDestination.hh 66241 2014-05-16 00:00:42Z adotti $
//
//
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4ofstreamDestination.hh
// Sends G4cout and G4cerr to a file.
//
// To send G4cout to a file instantiate an instance of
//     G4CoutToFile
// To send G4cerr to a file use instead
//     G4CerrToFile
//
// Example:
//  #include "G4ofstreamDestination.hh"
//  #include "G4ios.hh"
//  int main(int,char**) {
//      G4UImanager::GetUIpointer();
//      G4iosInitialization();
//      G4CoutToFile myout("MyOut_out.log"); // G4cout to file
//      G4CerrToFile myout("MyOut_err.log"); // G4cerr to file
// ---------------------------------------------------------------
#ifndef G4OFSTREAMDESTINATION_HH
#define G4OFSTREAMDESTINATION_HH

#include "G4coutDestination.hh"
#include <fstream>


class G4ofstreamDestinationBase : public G4coutDestination
{
public:
    G4ofstreamDestinationBase(const G4String& defaultName , G4bool append);
    virtual ~G4ofstreamDestinationBase();
    
    virtual G4int ReceiveG4cout( const G4String& ) = 0;
    virtual G4int ReceiveG4cerr( const G4String& ) = 0;
    
    void SetFileName( const G4String& name , G4bool append = true );
    void Close();
    void Open();
private:
    G4String fileName;
    G4bool appendFlag;
protected:
    std::ofstream g4file;
};

class G4CoutToFile : public G4ofstreamDestinationBase
{
public:
    G4CoutToFile(const G4String& deafaultName = "Geant4_cout.txt", G4bool append = true);
    virtual G4int ReceiveG4cout( const G4String& msg );
    virtual G4int ReceiveG4cerr( const G4String& ) { return 0;}
};

class G4CerrToFile : public G4ofstreamDestinationBase
{
public:
    G4CerrToFile(const G4String& deafaultName = "Geant4_cerr.txt", G4bool append = true);
    virtual G4int ReceiveG4cout( const G4String& ) { return 0; }
    virtual G4int ReceiveG4cerr( const G4String& msg );
    
};
#endif //G4OFSTREAMDESTINATION_HH
