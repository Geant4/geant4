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
// $Id: G4ITType.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4ITTYPE_HH
#define G4ITTYPE_HH 1

#include <cstddef>

/**
  * Tag the G4IT
  * Should be automatically setup by G4IT
  * using : ITDef(MyIT) and ITImp(MyIT)
  */

struct G4ITType
{
private :
    friend G4ITType operator +(const G4ITType& left,const int& right);
    friend G4ITType operator -(const G4ITType& left,const int& right);
    int fValue;

public :

    static size_t size();

    G4ITType(const int d_ = 0) : fValue(d_) {;}
    G4ITType(const G4ITType & d_) : fValue(d_.fValue){;}
    G4ITType & operator=(const G4ITType & rhs);
    inline G4ITType & operator=(const int & rhs) { fValue = rhs; return *this;}
    inline operator int & () { return fValue; }
    inline operator const int & () const { return fValue; }
    inline bool operator==(const G4ITType & rhs) const { return fValue == rhs.fValue; }
    inline bool operator==(const int & rhs) const { return fValue == rhs; }
    inline bool operator<(const G4ITType & rhs) const { return fValue < rhs.fValue; }
    inline void operator++() { fValue++; }
};

inline G4ITType operator +(const G4ITType& left,const int& right) {
    G4ITType output( left.fValue + right );
    return output;
}

inline G4ITType operator -(const G4ITType& left,const int& right) {
    G4ITType output( left.fValue - right );
    return output;
}

class G4ITTypeManager
{
private:
    static G4ITTypeManager* fgInstance ;
    G4ITType fLastType;
    G4ITTypeManager();
    virtual ~G4ITTypeManager();

public :
    G4ITType NewType() ;
    size_t size() const;
    static G4ITTypeManager* Instance();
    static void DeleteInstance();
};

#define ITDef(T)\
public:\
static G4ITType fType;\
static const G4ITType ITType()\
{\
    return fType;\
}\
const G4ITType GetITType() const\
{\
    return fType;\
}\
virtual G4bool equal(const G4IT &right) const \
{\
    const T& right_mol = (const T&)right ;\
    return (this->operator==(right_mol));\
}\
virtual G4bool diff(const G4IT &right) const\
{\
    const T& right_mol = (const T&)right ;\
    return (this->operator<(right_mol));\
}

#define ITImp(T) \
G4ITType T::fType = G4ITTypeManager::Instance()->NewType();

#endif // G4ITTYPE_HH
