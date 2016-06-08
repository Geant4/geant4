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
// $Id: defs.h,v 1.8.4.1 2001/06/28 19:09:57 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4RWBoundsErr, G4RWGeneralException
//
//  Class description:
//
//  STL wrapper classes for Errors and Exceptions utilities.
//  It implements Rogue Wave RWBoundsErr and RWGeneralException
//  signatures but intrinsically decoupled from Rogue Wave.

//---------------------------------------------------------------

#ifndef __defs
#define __defs

#include <stdio.h>
#include <string>
#include "G4Types.hh"

#define G4RWDEFAULT_CAPACITY (64)

#ifndef FALSE
  #define FALSE 0
#endif
#ifndef TRUE
  #define TRUE 1
#endif

#define CLHEP_MAX_MIN_DEFINED
#include <CLHEP/config/TemplateFunctions.h>

class G4RWBoundsErr
{
public:

  G4RWBoundsErr(const char* s,int b=0,int v=0)
    {
      char btmp[80],vtmp[80];
      sprintf(btmp,"%d",b);
      sprintf(vtmp,"%d",v);
      str=new char[strlen(s)+8+strlen(btmp)+7+strlen(vtmp)+1];
      strcpy(str,s);
      strcat(str,": bound:");
      strcat(str,btmp);
      strcat(str," value:");
      strcat(str,vtmp);
    }

  G4RWBoundsErr(const G4RWBoundsErr&e)
    {
      delete [] str;
      str=new char[strlen(e.str)+1];
      strcpy(str,e.str);
    }

  ~G4RWBoundsErr()
    {
      delete [] str;
    }

  const char * why() const
    {
      return str;
    }

private:

  char* str;

};
      
class G4RWGeneralException
{
public:

  G4RWGeneralException(const char* s)
    {
      str=new char[strlen(s)+1];
      strcpy(str,s);
    }

  G4RWGeneralException(const G4RWGeneralException&e)
    {
      delete [] str;
      str=new char[strlen(e.str)+1];
      strcpy(str,e.str);
    }

  ~G4RWGeneralException()
    {
      delete [] str;
    }

  const char * why() const
    {
      return str;
    }

private:

  char* str;

};
  
#ifndef G4NO_STD_EXCEPTIONS
  #define G4RWTHROW(a) throw a
#else
  #define G4RWTHROW(a) abort()
#endif

#endif
