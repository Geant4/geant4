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
// $Id: G4VRMLNetConfig.hh,v 1.4 2001-07-11 10:09:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRMLNetConfig.hh
// Satoshi Tanaka

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML_NET_CONFIG_HH
#define G4VRML_NET_CONFIG_HH

const int   FR_VRML_DEFAULT_PORT     = 40801              ;
const char  FR_VRML_PORT_ENV     []  = "G4VRML_PORT"      ;
const char  FR_VRML_HOST_NAME_ENV[]  = "G4VRML_HOST_NAME" ;

#endif // G4VRML_NET_CONFIG_HH
#endif // G4VIS_BUILD_VRML_DRIVER
