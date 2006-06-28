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
// File name:     RadmonDetectorFlatVolumeComponent.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeComponent.hh,v 1.3 2006-06-28 13:46:46 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Component to create a box
//

#ifndef   RADMONDETECTORFLATVOLUMECOMPONENT_HH
 #define  RADMONDETECTORFLATVOLUMECOMPONENT_HH
 
 // Include files
 #include "globals.hh"
 
 // Forward declarations
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4Box;
 class G4VisAttributes;
 
 class RadmonDetectorFlatVolumeComponent
 {
  public:
   inline                                       RadmonDetectorFlatVolumeComponent(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonDetectorFlatVolumeComponent();
   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeComponent();
                                                RadmonDetectorFlatVolumeComponent(const RadmonDetectorFlatVolumeComponent & copy);
   RadmonDetectorFlatVolumeComponent &          operator=(const RadmonDetectorFlatVolumeComponent & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   G4Box *                                      box;
   G4VisAttributes *                            visAttributes;
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeComponent.icc"
#endif /* RADMONDETECTORFLATVOLUMECOMPONENT_HH */
