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
// $Id: boxtoincpp.cc,v 1.3 2001-07-11 09:59:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ios.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"

const G4int norept=1000000;

void main()
{
    G4int i;
    G4double snxt,total=0;
    G4ThreeVector endpt,startpt(-500,0,0),dir(1,0,0);

    G4Box testbox(100,200,400);

    for (i=1;i<=norept;i++)
        {
            snxt=testbox.DistanceToIn(startpt,dir);
            endpt=G4ThreeVector(startpt.x()+snxt*dir.x(),
                                startpt.y()+snxt*dir.y(),
                                startpt.z()+snxt*dir.z());
            total=total+endpt.x()+endpt.y()+endpt.z();
        }    
    if ((!total)&&(!endpt.x())&&(!endpt.y())&&(!endpt.z()))
        {
            G4cout << "total=0";
        }
    return;
}
