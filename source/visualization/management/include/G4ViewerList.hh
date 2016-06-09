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
// $Id: G4ViewerList.hh,v 1.7 2003/06/16 17:14:12 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// 
// John Allison  May 1996

#ifndef G4VIEWERLIST_HH
#define G4VIEWERLIST_HH

#include <vector>
#include "G4VViewer.hh"

class G4ViewerList: public std::vector<G4VViewer*> {
public:
  void remove(G4VViewer*);
};

typedef std::vector<G4VViewer*>::iterator G4ViewerListIterator;
typedef std::vector<G4VViewer*>::const_iterator G4ViewerListConstIterator;

#endif
