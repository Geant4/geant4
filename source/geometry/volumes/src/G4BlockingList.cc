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
// $Id: G4BlockingList.cc,v 1.5 2002-04-19 08:22:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4BlockingList Implementation
//

#include "G4BlockingList.hh"

G4BlockingList::G4BlockingList(G4int maxDefault,G4int stride)
  : fStride(stride), fBlockingList(maxDefault,0)
{
}

// Do nothing destructor
G4BlockingList::~G4BlockingList()
{
}

// Clear List and reset tag
// Fix: Out of line for HP-CC
void G4BlockingList::FullyReset()
{
	fBlockTagNo=1;
	for (G4int i=fBlockingList.size()-1;i>=0;i--)
		{
			fBlockingList[i]=0;
		}	
	
}
