//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.
// Top level shell class

#ifndef G4VULKANXMTOPLEVELSHELL_HH
#define G4VULKANXMTOPLEVELSHELL_HH

#include "G4VulkanXmVWidgetShell.hh"

class G4VulkanXmVWidgetContainer;

class G4VulkanXmTopLevelShell : public G4VulkanXmVWidgetShell
{
  public:

    G4VulkanXmTopLevelShell(G4VulkanXmViewer*, char*);  // constructor
    virtual ~G4VulkanXmTopLevelShell();  // destructor

    void AddChild(G4VulkanXmVWidgetContainer*);
    void Realize();

    Widget* GetPointerToWidget();
    char* GetName();

  private:

    G4VulkanXmTopLevelShell(const G4VulkanXmTopLevelShell&);
    G4VulkanXmTopLevelShell& operator=(const G4VulkanXmTopLevelShell&);
    char* name;
    Widget toplevel;
    Widget top_box;
    Widget frame;
};

#endif
