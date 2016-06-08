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
// $Id: G4oCommon.hh,v 1.3.4.1 2001/06/28 19:06:32 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
/*
   Included by G4o.h.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*G4*/
#include <G4UImanager.hh>
#include <G4UIcommand.hh>
#include <G4UIparameter.hh>
#include <G4UIcommandTree.hh>

/*Co*/
#include <CMemory.h>
#include <CPrinter.h>
#include <CText.h>
#include <CFile.h>
#include <OType.h>
#include <OShell.h>
#include <OProcess.h>


static void         SetTypes                      ();
static void         ExecuteScript                 (char*);
static int          Execute_G4                    (int,char**,OProcess);
static int          DoMethod                      (OIdentifier,char*,void*,int,char**,void*,int*);
static G4UIcommand* GetCommandIdentifier          (char*);
static void         GetCommands                   (G4UIcommandTree*,int*,char***);
static int          ProduceODB                    (G4UIcommand*,char*);

