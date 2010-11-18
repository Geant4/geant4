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
#include <TCanvas.h>
#include <TH1.h>
#include <TAxis.h>
#include <TList.h>
#include "CexmcHistoWidget.hh"


void  CexmcHistoWidget::contextMenuEvent( QContextMenuEvent *  event )
{
    /* TRootContextMenu responsiveness is very ugly and unstable in Geant4
     * environment, so we have to replace it with something more useful
     * (e.g. unzooming histogram axes).
     * We can switcn back to context menu as soon as ROOT will implement Qt
     * backend for its context menu.
     * TQtWidget utilizes QContextMenuEvent::Other reason to draw
     * TRootContextMenu object. */
    if ( event->reason() == QContextMenuEvent::Other )
    {
        TCanvas *  canvas( Canvas() );
        TObject *  object( NULL );
        TListIter  curPrimitive( canvas->GetListOfPrimitives() );
        while ( ( object = curPrimitive.Next() ) )
        {
            if ( ! object->InheritsFrom( "TH1" ) )
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
        canvas->Update();
    }
}

#endif

