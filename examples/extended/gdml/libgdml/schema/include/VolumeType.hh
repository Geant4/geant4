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
// $Id: VolumeType.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef VOLUMETYPE_H
#define VOLUMETYPE_H 1

#include "IdentifiableVolumeType.hh"
#include "ReferenceType.hh"
#include "SinglePlacementType.hh"
#include "ContentGroup.hh"

class VolumeType : public IdentifiableVolumeType {
public:
    
  class materialref : public SAXObject, public ReferenceType {
  public:
    materialref() {
    }
    virtual ~materialref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class solidref : public SAXObject, public ReferenceType {
  public:
    solidref() {
    }
    virtual ~solidref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
public:
  VolumeType() {
  }
  ~VolumeType() {
  }
  
  const ContentSequence* get_content() const {
    return &m_sequence;
  }

  void add_content( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = { tag, so };
    m_sequence.add_content( ci );
  }
private:
  ContentSequence m_sequence;
};

#endif // VOLUMETYPE_H
