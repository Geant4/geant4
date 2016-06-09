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
// $Id: G4VEvaporationFactory.hh,v 1.2 2005/06/04 13:26:42 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4VEvaporationFactory_hh
#define G4VEvaporationFactory_hh


#include "G4VEvaporationChannel.hh"
#include <vector>

class G4VEvaporationFactory
{
public:
  G4VEvaporationFactory() : _channel(0) {};
  virtual ~G4VEvaporationFactory();

private:
  G4VEvaporationFactory(const G4VEvaporationFactory & ) {};
  const G4VEvaporationFactory & operator=(const G4VEvaporationFactory & val);
  G4bool operator==(const G4VEvaporationFactory & val) const;
  G4bool operator!=(const G4VEvaporationFactory & val) const;

public:
  
  std::vector<G4VEvaporationChannel*> * GetChannel();

protected:
  virtual std::vector<G4VEvaporationChannel*> * CreateChannel() = 0;



private:
  std::vector<G4VEvaporationChannel*> * _channel;

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
