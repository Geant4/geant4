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
// $Id: G4VTrajectoryPoint.hh,v 1.5 2001-07-11 10:08:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class description
//  The base class for the trajectory class

#ifndef G4VTrajectoryPoint_h
#define G4VTrajectoryPoint_h 1

class G4VTrajectoryPoint
{
   public:

   G4VTrajectoryPoint() {;}
   virtual ~G4VTrajectoryPoint() {;}

   inline int operator==(const G4VTrajectoryPoint& right) const
   { return (this==&right); }
};


#endif

