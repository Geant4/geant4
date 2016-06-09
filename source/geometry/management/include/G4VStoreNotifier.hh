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
// $Id: G4VStoreNotifier.hh,v 1.1 2004/09/02 07:48:58 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// class G4VStoreNotifier
//
// Class description:
//
// Simple abstract class allowing for implementation of user notifiers
// to be activated at registration/deregistration of objects in the
// volume, solid and region stores.

// Author:
// 01.09.04 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4VSTORENOTIFIER_HH
#define G4VSTORENOTIFIER_HH

class G4VStoreNotifier
{
  public:  // with description

    G4VStoreNotifier();
    virtual ~G4VStoreNotifier();
      // Constructor and destructor.

    virtual void NotifyRegistration() = 0;
      // Notification of object registration.
    virtual void NotifyDeRegistration() = 0;
      // Notification of object deregistration.
};

#endif
