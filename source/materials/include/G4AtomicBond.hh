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
//
//
//---------------------------------------------------------------
//
// G4AtomicBond
//
// Class Description:
//


#ifndef G4AtomicBond_H
#define G4AtomicBond_H 1

#include "globals.hh"
#include "G4Element.hh"

class G4AtomicBond
{
public:
    enum theBondType {Ionic=0, Covalent=1, Metallic=2, NA=-1 };

    G4AtomicBond(theBondType aType,
                 G4Element* firstAtomKind,
                 G4int firstAtomNumber,
                 G4Element* secondAtomKind,
                 G4int secondAtomNumber);
    
    virtual ~G4AtomicBond();
    
    //
    // The indexes of the two atoms connected by the bond are stored
    // theFirstAtomKind is the index in the vector of G4Elements in the array
    // and theFirstAtomNumber is the index in the vector of G4ThreeVector in
    // the G4CrystalLattice. For the second atom theSecondAtomKind and
    // theSecondAtomNumber are used.
    //
private:
    G4Element* theFirstAtomKind;
    G4int theFirstAtomNumber;

    G4Element* theSecondAtomKind;
    G4int theSecondAtomNumber;

public:
    inline G4Element* GetFirstAtomKind() const {return theFirstAtomKind;};
    void SetFirstAtomKind(G4Element* aElement)  {theFirstAtomKind=aElement;};

    inline G4int GetFirstAtomNumber() const {return theFirstAtomNumber;};
    void SetFirstAtomNumber(G4int aInt)  {theFirstAtomNumber=aInt;};

    inline G4Element* GetSecondAtomKind() const {return theSecondAtomKind;};
    void SetSecondAtomKind(G4Element* aElement)  {theSecondAtomKind=aElement;};
    
    inline G4int GetSecondAtomNumber() const {return theSecondAtomNumber;};
    void SetSecondAtomNumber(G4int aInt)  {theSecondAtomNumber=aInt;};

    //
    // theType stores which kind of bond is the G4AtomicBond object
    // there are three types: Ionic, Covalent and Metallic
    // and they are enumerate in the theBondType enum.
private:
    theBondType theType;

public:
    inline theBondType GetType() const {return theType;};
    void SetType(theBondType aType)  {theType=aType;};

    //
    // theAromaticity stores if the bond is aromatic.
    //
private:
    G4bool theAromaticity;

public:
    inline G4bool GetAromaticity() const {return theAromaticity;};
    void SetAromaticity(G4bool aBool)  {theAromaticity=aBool;};

    
};

#endif

