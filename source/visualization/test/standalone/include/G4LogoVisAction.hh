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
// $Id: G4LogoVisAction.hh,v 1.2 2005-03-16 17:25:15 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4LOGOVISACTION_HH
#define G4LOGOVISACTION_HH

#include "G4VUserVisAction.hh"

class G4VisAttributes;
class G4Polyhedron;

class G4LogoVisAction: public G4VUserVisAction {
public:
  G4LogoVisAction();
  ~G4LogoVisAction();
  void Draw();
private:
  G4VisAttributes* fpVisAtts;
  G4Polyhedron *fpG, *fp4;
};

#endif

