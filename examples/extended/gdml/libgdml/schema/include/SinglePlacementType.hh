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
// $Id: SinglePlacementType.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SINGLEPLACEMENTTYPE_H
#define SINGLEPLACEMENTTYPE_H 1

#include "ContentGroup.hh"
#include "ReferenceType.hh"

class SinglePlacementType {
public:
    
  class volumeref : public SAXObject, public ReferenceType {
  public:
    volumeref() {
    }
    virtual ~volumeref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
//  We need to resolve the problem of conlifting non-terminals in grammar
//  This is typical use-case here, the following two elements clash with
//  BooleanSolidType elements
//   class positionref : public SAXObject, public ReferenceType {
//   public:
//     positionref() {
//     }
//     virtual ~positionref() {
//     }
//     virtual SAXObject::Type type() {
//       return SAXObject::element;
//     }
//   };
//   
//   class rotationref : public SAXObject, public ReferenceType {
//   public:
//     rotationref() {
//     }
//     virtual ~rotationref() {
//     }
//     virtual SAXObject::Type type() {
//       return SAXObject::element;
//     }
//   };
  
public:
  SinglePlacementType() {
  }
  ~SinglePlacementType() {
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

#endif // SINGLEPLACEMENTTYPE_H
