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
// $Id: G4VNotifier.hh,v 1.1 2006-05-03 09:59:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VNotifier
//
// Class description:
//
// Simple abstract class allowing for implementation of user notifiers
// to be activated for example at registration/deregistration of objects
// in stores.

// Author:
// 01.09.04 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4VNOTIFIER_HH
#define G4VNOTIFIER_HH

class G4VNotifier
{
  public:  // with description

    G4VNotifier();
    virtual ~G4VNotifier();
      // Constructor and destructor.

    virtual void NotifyRegistration() = 0;
      // Notification of object registration.
    virtual void NotifyDeRegistration() = 0;
      // Notification of object deregistration.
};

#endif
