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
// File name:     RadmonRunActionObserver.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonRunActionObserver.hh,v 1.2 2006-06-28 13:57:27 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class for the Run Action
//

#ifndef   RADMONRUNACTIONOBSERVER_HH
 #define  RADMONRUNACTIONOBSERVER_HH

 class G4Run;
 
 class RadmonRunActionObserver
 {
  public:
   inline                                       RadmonRunActionObserver();
   inline virtual                              ~RadmonRunActionObserver();
   
   virtual void                                 OnBeginOfRun(const G4Run * run) = 0;
   virtual void                                 OnEndOfRun(const G4Run * run) = 0;

  private:
                                                RadmonRunActionObserver(const RadmonRunActionObserver & copy);
   RadmonRunActionObserver &                    operator=(const RadmonRunActionObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonRunActionObserver.icc"
#endif /* RADMONRUNACTIONOBSERVER_HH */
