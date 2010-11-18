/*
 * =============================================================================
 *
 *       Filename:  CexmcHistoWidget.hh
 *
 *    Description:  histogram widget without context menu
 *                  (derived from TQtWidget)
 *
 *        Version:  1.0
 *        Created:  15.03.2010 18:59:09
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HISTO_WIDGET_HH
#define CEXMC_HISTO_WIDGET_HH

#ifdef CEXMC_USE_ROOTQT

#include <TQtWidget.h>

class  QContextMenuEvent;


class  CexmcHistoWidget : public TQtWidget
{
    public:
        void  contextMenuEvent( QContextMenuEvent *  event );
};

#endif

#endif

