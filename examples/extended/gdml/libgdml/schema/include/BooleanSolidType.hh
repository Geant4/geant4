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
// $Id: BooleanSolidType.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef BOOLEANSOLIDTYPE_H
#define BOOLEANSOLIDTYPE_H 1

#include "SAXObject.hh"

#include "SolidType.hh"
#include "ContentGroup.hh"
#include "ReferenceType.hh"
#include "positionType.hh"
#include "rotationType.hh"

class BooleanSolidType : public SolidType {
public:
  class first : public SAXObject, public ReferenceType {
  public:
    first() {
    }
    virtual ~first() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class second : public SAXObject, public ReferenceType {
  public:
    second() {
    }
    virtual ~second() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class position : public SAXObject, public positionType {
  public:
    position() {
    }
    virtual ~position() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class positionref : public SAXObject, public ReferenceType {
  public:
    positionref() {
    }
    virtual ~positionref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class rotation : public SAXObject, public rotationType {
  public:
    rotation() {
    }
    virtual ~rotation() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class rotationref : public SAXObject, public ReferenceType {
  public:
    rotationref() {
    }
    virtual ~rotationref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
public:
  BooleanSolidType() {
  }
  
  ~BooleanSolidType() {
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

#endif // BOOLEANSOLIDTYPE_H
