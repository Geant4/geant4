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
 * File:   G4ArrayOps.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 28, 2012, 16:08
 */

// All the arrays that use the functions in this namespace MUST have a
// [] (bracket) operator for element-wise access

#ifndef G4ARRAYOPS_HH
#define	G4ARRAYOPS_HH

#include <vector>

#include "globals.hh"

/** G4ArrayOps is a namespace that provides template functions for
 *  performing basic arithmatic operations on any data type that
 *  is accessed with the [] operator.
 */
namespace G4ArrayOps
{
    /** Set's all the values in an array to a constant */
    template< class T >
    void Set( G4int Elements,
               T* To,
               T Value )
    {
        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = Value;
        }
    }

    /** Copy values from one array to another */
    template< class T >
    void Copy( G4int Elements,
               T* To,
               T* From )
    {
        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = From[position];
        }
    }

    /** Add two arrays together. If the second array is NULL then the
     *  'To' array is used as if the function were the += operator.
     */
    template< class T >
    void Add( G4int Elements,
              T* To,
              T* A1,
              T* A2 = NULL )
    {
        if(A2 == NULL)
        {
            A2 = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = A2[position] + A1[position];
        }
    }

    /** Add a constant to an array. If the second array is NULL then the
     *  'To' array is used as if the function were the += operator.
     */
    template< class T >
    void Add( G4int Elements,
              T* To,
              T A1,
              T* A2 = NULL )
    {
        if(A2 == NULL)
        {
            A2 = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = A1 + A2[position];
        }
    }

    /** Subtract an array from another. If the second array is NULL then the
     *  'To' array is used as if the function were the -= operator.
     */
    template< class T >
    void Subtract( G4int Elements,
                   T* To,
                   T* Minuend,
                   T* Subtrahend = NULL )
    {
        if(Subtrahend == NULL)
        {
            Subtrahend = Minuend;
            Minuend = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = Minuend[position] - Subtrahend[position];
        }
    }

    /** Multiply two arrays together. If the second array is NULL then the
     *  'To' array is used as if the function were the *= operator.
     */
    template< class T >
    void Multiply( G4int Elements,
                   T* To,
                   T* M1,
                   T* M2 = NULL )
    {
        if(M2 == NULL)
        {
            M2 = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = M2[position] * M1[position];
        }
    }

    /** Multiply an array by a constant. If the second array is NULL then the
     *  'To' array is used as if the function were the *= operator.
     */
    template< class T >
    void Multiply( G4int Elements,
                   T* To,
                   T M1,
                   T* M2 = NULL )
    {
        if(M2 == NULL)
        {
            M2 = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = M2[position] * M1;
        }
    }

    /** Divide an array by another. If the second array is NULL then the
     *  'To' array is used as if the function were the /= operator.
     */
    template< class T >
    void Divide( G4int Elements,
                 T* To,
                 T* Numerator,
                 T* Denominator = NULL )
    {
        if(Denominator == NULL)
        {
            Denominator = Numerator;
            Numerator = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = Numerator[position] / Denominator[position];
        }
    }

    /** Divide a constant by an array. If the second array is NULL then the
     *  'To' array is used as if the function were the /= operator.
     */
    template< class T >
    void Divide( G4int Elements,
                 T* To,
                 T Numerator,
                 T* Denominator = NULL )
    {
        if(Denominator == NULL)
        {
            Denominator = To;
        }

        for(G4int position = 0; position < Elements; position++)
        {
            To[position] = Numerator / Denominator[position];
        }
    }

    template< class T >
    void DeleteVectorOfPointers( std::vector< T >& Vector )
    {
    	for(unsigned int i = 0; i < Vector.size(); i++)
    	{
			delete Vector[i];
		}

    	delete &Vector;
    }
}

#endif  /** G4ARRAYOPS_HH */

