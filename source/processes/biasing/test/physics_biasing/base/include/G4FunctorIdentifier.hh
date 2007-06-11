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
// $Id: G4FunctorIdentifier.hh,v 1.2 2007-06-11 19:25:47 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation.
//
#ifndef G4FUNCTORIDENTIFIER_HH
#define G4FUNCTORIDENTIFIER_HH

#include "G4String.hh"

class G4FunctorIdentifier {

public:
  typedef unsigned long Key;
  typedef G4String StringType;

  G4FunctorIdentifier(const StringType& name = "null",
		      Key parentID=0)
  {
    fParentID = parentID;
    fUniqueID = NextKey();
    fName = name;
  }

  Key UniqueID() const {return fUniqueID;}
  Key ParentID() const {return fParentID;}
  const StringType& Name() const {return fName;}

  bool operator == (const G4FunctorIdentifier& other) const
  {
    return ((fName == other.fName) &&
            (fUniqueID == other.fUniqueID) &&
            (fParentID == other.fParentID));
  }
private:

  Key NextKey() const {
    static Key nKey = 0;
    return ++nKey;
  }

  StringType fName;
  Key fUniqueID;
  Key fParentID;

};

#endif

// jane fixme: to do
// 1) What to use for name - G4String, std::string or what ? speed issues with G4String apparently

