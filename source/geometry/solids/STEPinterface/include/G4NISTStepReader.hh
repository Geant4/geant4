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
// $Id: G4NISTStepReader.hh,v 1.5 2001-07-11 10:00:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4NISTStepReader
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4NISTSTEPFILEREADER_HH
#define G4NISTSTEPFILEREADER_HH

#include <instmgr.h>
#include "G4StepFileReader.hh"

///////////////////////////////////////////////////////////////////////////////
// Function defined as a stub (necessary to use the scl)
// The purpose of this function is to allow the DisplayNode object to delete 
// an object that it knows nothing about.  It was made generic so that the scl
// could be used with any display toolkit.
//
// This function is called by the DisplayNode object
// This function needs to be defined outside the SCL libraries.  It needs to do
// two things:
// 1) unmap the StepEntityEditor window if it is mapped.
// 2) delete the StepEntityEditor window
// To see an example of this function used with the Data Probe look in
// ../clprobe-ui/StepEntEditor.cc  Look at DeleteSEE() and ~StepEntityEditor().
///////////////////////////////////////////////////////////////////////////////

// this function illustrates a good way to Generate and assign file identifiers
/*
void AssignFileId (STEPentity *se, InstMgr& instance_list)
{
    int fId = instance_list.MaxFileId() + 1;
    se->STEPfile_id = (fId > 0) ? fId : 1; 
}
*/

// define this to be the name of the display window object for 
// STEP entity instance editing or define your own.
// This is only needed as there's a link to these from the toolkit

class StepEntityEditor
{
  public:
    StepEntityEditor() {;}
    ~StepEntityEditor() {;}
};

extern void AssignFileId (STEPentity *se, InstMgr& instance_list);
extern STEPentity *GetEntity (STEPnode *node, InstMgr *im);

// This needs to be defined for the STEPfile reader
extern void SchemaInit (Registry &);

class G4NISTStepReader : public G4StepFileReader
{
  public:

    void ReadSTEPFile(const char*);
    void SaveSTEPFile();
    void UpdateSTEPFile();
    InstMgr GetInstanceManager() const { return InstanceList; }

  private:

    InstMgr InstanceList;
};

#endif
