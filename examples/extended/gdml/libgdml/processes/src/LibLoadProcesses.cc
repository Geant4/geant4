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
// $Id: LibLoadProcesses.cc,v 1.2 2002-06-03 12:09:32 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "SAXComponentFactory.hh"

extern "C" {

  void GDMLProcessLibLoad() {
    LOAD_COMPONENT(defineProcess)
    LOAD_COMPONENT(constantProcess)
    LOAD_COMPONENT(quantityProcess)
    LOAD_COMPONENT(expressionProcess)
    LOAD_COMPONENT(positionProcess)
    LOAD_COMPONENT(rotationProcess)
    LOAD_COMPONENT(atomProcess)
    LOAD_COMPONENT(DProcess)
    LOAD_COMPONENT(DrefProcess)
    LOAD_COMPONENT(TProcess)
    LOAD_COMPONENT(TrefProcess)
    LOAD_COMPONENT(PProcess)
    LOAD_COMPONENT(PrefProcess)
    LOAD_COMPONENT(isotopeProcess)
    LOAD_COMPONENT(elementProcess)
    LOAD_COMPONENT(materialProcess)
    LOAD_COMPONENT(fractionProcess)
    LOAD_COMPONENT(compositeProcess)    
    LOAD_COMPONENT(firstProcess)    
    LOAD_COMPONENT(secondProcess)    
    LOAD_COMPONENT(positionrefProcess)    
    LOAD_COMPONENT(rotationrefProcess)    
    LOAD_COMPONENT(unionProcess)    
    LOAD_COMPONENT(subtractionProcess)    
    LOAD_COMPONENT(intersectionProcess)
    LOAD_COMPONENT(boxProcess)    
    LOAD_COMPONENT(sphereProcess)    
    LOAD_COMPONENT(tubeProcess)    
    LOAD_COMPONENT(coneProcess)    
    LOAD_COMPONENT(paraProcess)    
    LOAD_COMPONENT(trdProcess)    
    LOAD_COMPONENT(trapProcess)
    LOAD_COMPONENT(volumeProcess)
    LOAD_COMPONENT(assemblyProcess)
    LOAD_COMPONENT(childProcess)
    LOAD_COMPONENT(materialrefProcess)    
    LOAD_COMPONENT(solidrefProcess)    
    LOAD_COMPONENT(volumerefProcess)
    LOAD_COMPONENT(worldProcess)    
    LOAD_COMPONENT(setupProcess)    
  }

};
