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

#ifndef G4XNNTotal_h
#define G4XNNTotal_h

#include "globals.hh"
#include "G4CrossSectionPatch.hh"
#include "G4CrossSectionVector.hh"

class G4KineticTrack;

class G4XNNTotal : public G4CrossSectionPatch
{

public:

  G4XNNTotal();

  virtual ~G4XNNTotal();

  G4bool operator==(const G4XNNTotal &right) const;
  G4bool operator!=(const G4XNNTotal &right) const;
  virtual const G4CrossSectionVector* GetComponents() const { return components; } 
 
  virtual G4String Name() const;


protected:


private:  

  G4XNNTotal(const G4XNNTotal &right);
  const G4XNNTotal& operator=(const G4XNNTotal &right);

  G4CrossSectionVector* components;

};

#endif
















