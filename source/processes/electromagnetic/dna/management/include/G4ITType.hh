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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef G4ITTYPE_HH
#define G4ITTYPE_HH 1

#include <cstddef>
#include "G4Types.hh"

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
    inline G4bool operator==(const G4ITType & rhs) const { return fValue == rhs.fValue; }
    inline G4bool operator==(const int & rhs) const { return fValue == rhs; }
    inline G4bool operator<(const G4ITType & rhs) const { return fValue < rhs.fValue; }
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
    static /*G4ThreadLocal*/ G4ITTypeManager* fgInstance ;
    static G4ThreadLocal G4ITTypeManager* fgInstance_local ;
    G4ITType fLastType;
    G4ITTypeManager();
    virtual ~G4ITTypeManager();

    size_t fRessource;

public :
    G4ITType NewType() ;
    size_t size() const;
    static G4ITTypeManager* Instance();
    static void DeleteInstance();

    void ReserveRessource();
    void ReleaseRessource();
};

#define ITDef(T)\
public:\
inline static G4ITType fType= G4ITTypeManager::Instance()->NewType();\
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

#endif // G4ITTYPE_HH
