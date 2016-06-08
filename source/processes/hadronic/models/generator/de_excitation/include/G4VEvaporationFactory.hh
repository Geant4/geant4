//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VEvaporationFactory.hh,v 1.3 2002/12/12 19:17:14 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4VEvaporationFactory_hh
#define G4VEvaporationFactory_hh


#include "G4VEvaporationChannel.hh"
#include "g4std/vector"

class G4VEvaporationFactory
{
public:
  G4VEvaporationFactory() : _channel(0) {};
  virtual ~G4VEvaporationFactory();

private:
  G4VEvaporationFactory(const G4VEvaporationFactory & val) {};
  const G4VEvaporationFactory & operator=(const G4VEvaporationFactory & val);
  G4bool operator==(const G4VEvaporationFactory & val) const;
  G4bool operator!=(const G4VEvaporationFactory & val) const;

public:
  
  G4std::vector<G4VEvaporationChannel*> * GetChannel();

protected:
  virtual G4std::vector<G4VEvaporationChannel*> * CreateChannel() = 0;



private:
  G4std::vector<G4VEvaporationChannel*> * _channel;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };
  
};
#endif
