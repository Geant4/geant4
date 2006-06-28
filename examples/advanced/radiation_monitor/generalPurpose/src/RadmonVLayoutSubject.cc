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
// File name:     RadmonVLayoutSubject.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVLayoutSubject.cc,v 1.2 2006-06-28 13:52:41 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonVLayoutSubject.hh"

#include "RadmonVLayoutObserver.hh"


void                                            RadmonVLayoutSubject :: AttachObserver(RadmonVLayoutObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonVLayoutSubject :: DetachObserver(RadmonVLayoutObserver * observer)
{
 observersSet.erase(observer);
}





void                                            RadmonVLayoutSubject :: NotifyChange(void)
{
 ObserversSet::iterator i=observersSet.begin();
 ObserversSet::iterator end=observersSet.end();

 while (i!=end)
 {
  (*i)->OnLayoutChange();
  i++;
 }
}
