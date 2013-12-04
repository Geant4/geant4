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
/*
 * File:   G4TableTemplate.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on August 2, 2011, 10:21 AM
 */

#ifndef G4TABLETEMPLATE_HH
#define	G4TABLETEMPLATE_HH

#include <vector>

#include "globals.hh"

#include "G4FFGDefaultValues.hh"

/** G4TableTemplate is essentially a wrapper around a std::vector designed
 *  to work specifically with pointers.
 */
template <class T>
class G4TableTemplate
{
public:
    /** Default constructor */
    G4TableTemplate( void );
    /** Adds a container to the table */
    void G4AddContainer( T* NewContainer );
    /** Gets a pointer to the table */
    G4TableTemplate* G4GetTable( void );
    /** Retrieve a container from the table */
    T* G4GetContainer( unsigned int WhichContainer );
    /** Create a new blank container */
    T* G4GetNewContainer( void );
    /** Create a new container that is constructed with a G4int */
    T* G4GetNewContainer( G4int DefaultValue );
    /** Get the number of elements in the table */
    G4long G4GetNumberOfElements( void );

private:
    std::vector<T*> ContainerTable_;

public:
    ~G4TableTemplate( void );
};

template <class T>
G4TableTemplate<T>::
G4TableTemplate( void )
{
    // Nothing to be initialized
}

template <class T>
void G4TableTemplate<T>::
G4AddContainer( T* NewContainer )
{
    ContainerTable_.push_back(NewContainer);
}

template <class T>
G4TableTemplate<T>* G4TableTemplate<T>::
G4GetTable( void )
{
    return this;
}

template <class T>
T* G4TableTemplate<T>::
G4GetContainer( unsigned int WhichContainer )
{
    if(WhichContainer < ContainerTable_.size())
    {
        return ContainerTable_[WhichContainer];
    }

    return NULL;
}

template <class T>
T* G4TableTemplate<T>::
G4GetNewContainer( void )
{
    ContainerTable_.push_back(new T);

    return ContainerTable_.back();
}

template <class T>
T* G4TableTemplate<T>::
G4GetNewContainer( G4int DefaultValue )
{
    ContainerTable_.push_back(new T(DefaultValue));

    return ContainerTable_.back();
}

template <class T>
G4long G4TableTemplate<T>::
G4GetNumberOfElements( void )
{
    return ContainerTable_.size();
}

template <class T>
G4TableTemplate<T>::
~G4TableTemplate()
{
    for(unsigned int i = 0; i < ContainerTable_.size(); i++)
    {
        delete ContainerTable_[i];
    }
}

#endif	/* G4TABLETEMPLATE_HH */

