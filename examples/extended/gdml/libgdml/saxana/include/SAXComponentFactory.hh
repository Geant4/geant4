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
// $Id: SAXComponentFactory.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_COMPONENT_FACTORY_H
#define SAX_COMPONENT_FACTORY_H 1

#include "SAXComponentFactoryBase.hh"
#include "SAXComponentFactoryTable.hh"

template <class TYPE>
class SAXComponentFactory : virtual public SAXComponentFactoryBase
{
public:
  SAXComponentFactory( SAXComponentObject::EType type )
  : SAXComponentFactoryBase(), fType( type )
  {
    SAXComponentFactoryTable::GetInstance()->Register( this );
  }
  
  ~SAXComponentFactory()
  {
  }
  virtual SAXComponentObject* Build() const;
  virtual SAXComponentObject::EType Type() const;

private:
  SAXComponentObject::EType fType;
};

template <class TYPE>
SAXComponentObject* SAXComponentFactory<TYPE>::Build() const
{
  SAXComponentObject* obj = new TYPE();
  return obj;
}

template <class TYPE>
SAXComponentObject::EType SAXComponentFactory<TYPE>::Type() const
{
  return fType;
}

#define DECLARE_PROCESS_FACTORY(x) \
static const SAXComponentFactory<x> s_##x##Factory( SAXComponentObject::eProcess ); \
const SAXComponentFactoryBase& x##FactoryRef = s_##x##Factory;

#define DECLARE_SUBSCRIBER_FACTORY(x) \
static const SAXComponentFactory<x> s_##x##Factory( SAXComponentObject::eSubscriber ); \
const SAXComponentFactoryBase& x##FactoryRef = s_##x##Factory;

#define LOAD_COMPONENT(x)         extern const SAXComponentFactoryBase& x##FactoryRef; x##FactoryRef.Load();

#endif // SAX_COMPONENT_FACTORY_H

