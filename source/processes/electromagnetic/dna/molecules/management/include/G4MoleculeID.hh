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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef  G4MoleculeID_h
#define  G4MoleculeID_h 1

class G4MoleculeID
{
private :
    friend G4MoleculeID operator +(const G4MoleculeID& left,const int& right);
    friend G4MoleculeID operator -(const G4MoleculeID& left,const int& right);
    int fValue;

    static int fLastValue;

public :

    static int Last()
    {
        return fLastValue ;
    }

    static G4MoleculeID Create()
    {
        fLastValue ++;
        return G4MoleculeID(fLastValue);
    }

    static G4MoleculeID Initialize(int i)
    {
        fLastValue = i;
        return G4MoleculeID(i);
    }

    G4MoleculeID(const int d_) : fValue(d_){;}

    G4MoleculeID(){fValue=0;}
    G4MoleculeID(const G4MoleculeID & d_) : fValue(d_.fValue){;}
    inline G4MoleculeID & operator=(const G4MoleculeID & rhs) { this->fValue = rhs.fValue; return *this;}
    G4MoleculeID & operator=(const int & rhs) { this->fValue = rhs; return *this;}
    inline operator int & () { return fValue; }
    inline operator const int & () const { return fValue; }
    inline bool operator==(const G4MoleculeID & rhs) const { return fValue == rhs.fValue; }
    inline bool operator==(const int & rhs) const { return fValue == rhs; }
    inline bool operator<(const G4MoleculeID & rhs) const { return fValue < rhs.fValue; }
};

extern G4MoleculeID gStartCounter;

#endif
