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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4VITPROCESSOR_H
#define G4VITPROCESSOR_H

#include "globals.hh"
#include <stdlib.h>


template <typename PROC>
class G4VITProcessor
{
//--------
protected :
//--------
    // Protected member functions
    static void* start(void*);
    int start_thread(void*);

    typedef void (PROC::*Action)(void) ;
    void run(Action);

    G4VITProcessor(const G4VITProcessor& /*other*/);
    G4VITProcessor& operator=(const G4VITProcessor& /*other*/);
    template <typename U> friend void Join(G4VITProcessor<U>* processor);

//--------
public:
//--------
    // Constructors & Desctructors
    G4VITProcessor(G4int /*nbProc*/);
    virtual ~G4VITProcessor();

    // Public member methods

    virtual void Initialize(void* anObject = 0);
    // Becareful : use default constructor

    inline G4int        GetNbProc()            { return fNbProc ;  }
    inline bool         Occupied()
    {
        if(fNbProc == 0) return false ;
        G4bool output = fOccupied;

        if(output == false && fRunning == true)
        {
            pthread_join(fThreadID,NULL);
            fRunning = false;
        }
        return output;
    }
    inline G4int        GetID()                { return fID ;      }
    inline PROC*   GetNext()              { return fNext ;    }

    inline pthread_t& GetThreadID()
    {
        return fThreadID;
    }

//--------
protected:
//--------
    // Protected Attributes
    static int fNbProc ;
    bool fRunning ;
    //*********************************************************
    PROC* fNext ;
    int fID;
    // if fID == fNbProc => first object instanciated

    //*********************************************************
    pthread_t fThreadID;
    pthread_mutex_t fMutex;
    static pthread_attr_t fThreadAttribute;

    bool fInitialized;
    /*volatile*/ bool fOccupied ;

    //*********************************************************
    // Store method to apply a method
    // peut etre template ?
    struct ApplyMethod
    {
        PROC* fProcessor ;
        Action fAction ;
        G4bool fThreaded ;

        // Constructor & Destructor
        ApplyMethod(PROC* anObject, Action anAction)
        {
            fProcessor = anObject;
            fAction = anAction ;
            fThreaded = false;
        }

        ~ApplyMethod()
        {
            fProcessor = 0;
        }
    } ;

    //*********************************************************
    void printError(const char * msg, int status, const char* fileName, int lineNumber)
    {
        G4cerr << msg << " " << fileName << ":" << lineNumber <<
        "-" << strerror(status) << G4endl;
    }

//------
private:
//------

};

template< typename PROC>
void Join(G4VITProcessor<PROC>* processor)
{
    if(processor->GetNbProc() == 0) return ;
    G4int IDstart = processor->GetID() ;
    G4VITProcessor<PROC>* tmpProcessor = processor;
    void *status;

    do
    {
        if(tmpProcessor->fRunning)
        {
            int rc = pthread_join(tmpProcessor->GetThreadID(), &status);
            if (rc)
            {
                G4String orignException(__PRETTY_FUNCTION__);
                orignException += " ";
                orignException += __FILE__;
                orignException += " (";
                std::ostringstream os;
                os << __LINE__;
                orignException += os.str();
                orignException += ")";
                G4String exceptionCode ("ITProcessor001");
                G4ExceptionDescription exceptionDescription ;
                exceptionDescription << "return code from pthread_join() is :";
                exceptionDescription << rc ;
                G4Exception(orignException.data(),exceptionCode.data(),
                            FatalErrorInArgument,exceptionDescription);
            }
            tmpProcessor->fRunning = false;
        }
        tmpProcessor = tmpProcessor->GetNext();
    }
    while(tmpProcessor->GetID() != IDstart);
}

#ifdef TEMPLATE
#undef TEMPLATE
#endif

#define TEMPLATE template<typename PROC>
#define G4VITPROCESSOR G4VITProcessor<PROC>

#include "G4VITProcessor.icc"

#undef TEMPLATE
#undef G4VITPROCESSOR
#endif // G4VITPROCESSOR_H
