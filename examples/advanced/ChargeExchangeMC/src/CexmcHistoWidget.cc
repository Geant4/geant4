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
 * =============================================================================
 *
 *       Filename:  CexmcHistoWidget.cc
 *
 *    Description:  histogram widget without context menu
 *                  (derived from TQtWidget)
 *
 *        Version:  1.0
 *        Created:  15.03.2010 19:01:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifdef CEXMC_USE_ROOTQT

#include <QContextMenuEvent>
#include <TObject.h>
#include <TH1.h>
#include <TAxis.h>
#include <TList.h>
#include "CexmcHistoWidget.hh"


CexmcHistoWidget::CexmcHistoWidget()
{
    /* this is a workaround of the repaint bug in the ROOT Qt backend:
     * see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=17081#p73055 */
    fCanvas->SetFillColor( 10 );
}


void  CexmcHistoWidget::contextMenuEvent( QContextMenuEvent *  aEvent )
{
    /* TRootContextMenu responsiveness is very ugly and unstable in Geant4
     * environment, so we have to replace it with something more useful
     * (e.g. unzooming histogram axes).
     * We can switcn back to context menu as soon as ROOT will implement Qt
     * backend for its context menu.
     * TQtWidget utilizes QContextMenuEvent::Other reason to draw
     * TRootContextMenu object. */
    if ( aEvent->reason() == QContextMenuEvent::Other )
    {
        TObject *  object( NULL );
        TListIter  curPrimitive( fCanvas->GetListOfPrimitives() );
        while ( ( object = curPrimitive.Next() ) )
        {
            if ( ! object->InheritsFrom( TH1::Class() ) )
                continue;
            TH1 *    histo( static_cast< TH1 * >( object ) );
            TAxis *  axis( histo->GetXaxis() );
            if ( axis )
                axis->UnZoom();
            axis = histo->GetYaxis();
            if ( axis )
                axis->UnZoom();
            axis = histo->GetZaxis();
            if ( axis )
                axis->UnZoom();
        }
        fCanvas->Update();
    }
}

#endif

