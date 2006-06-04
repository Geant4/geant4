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
// $Id: ZBase.hh,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   ZBase.hh
//
//                                         2005 Q
// ====================================================================
#ifndef ZBASE_H
#define ZBASE_H

#include <string>

// ====================================================================
//
// class definition
//
// ====================================================================

class ZBase {
public:
  ZBase();
  virtual ~ZBase();

  void AMethod();
  virtual void VMethod(std::string message);

};

#endif
