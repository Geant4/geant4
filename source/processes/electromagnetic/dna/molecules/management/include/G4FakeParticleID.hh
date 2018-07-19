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


#ifndef  G4FakeParticleID_h
#define  G4FakeParticleID_h 1

#include "G4Types.hh" 

class G4FakeParticleID
{
private :
    friend G4FakeParticleID operator +(const G4FakeParticleID& left,const int& right);
    friend G4FakeParticleID operator -(const G4FakeParticleID& left,const int& right);
    int fValue;

    static G4ThreadLocal int fLastValue;

public :

    static int Last()
    {
        return fLastValue ;
    }

    static G4FakeParticleID Create()
    {
        fLastValue ++;
        return G4FakeParticleID(fLastValue);
    }

    static G4FakeParticleID Initialize(int i)
    {
        fLastValue = i;
        return G4FakeParticleID(i);
    }

    G4FakeParticleID(const int d_) : fValue(d_){;}

    G4FakeParticleID(){fValue=0;}
    G4FakeParticleID(const G4FakeParticleID & d_) : fValue(d_.fValue){;}
    inline G4FakeParticleID & operator=(const G4FakeParticleID & rhs) { this->fValue = rhs.fValue; return *this;}
    G4FakeParticleID & operator=(const int & rhs) { this->fValue = rhs; return *this;}
    inline operator int & () { return fValue; }
    inline operator const int & () const { return fValue; }
    inline bool operator==(const G4FakeParticleID & rhs) const { return fValue == rhs.fValue; }
    inline bool operator==(const int & rhs) const { return fValue == rhs; }
    inline bool operator<(const G4FakeParticleID & rhs) const { return fValue < rhs.fValue; }
};

extern G4FakeParticleID gStartCounter;

#endif
