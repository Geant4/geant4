/* file scrollmouse.c
 *
 * Event handler for wheel mouse
 *
 * functions: void xmAddMouseEventHandler(Widget w)
 *
 */
/*
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <Xm/Xm.h>
#include <Xm/ScrollBar.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/*******************************************************************/ 
/*                                                                 */ 
/* NAME:      mouseScroll                                          */ 
/*                                                                 */ 
/* FUNCTION:  do scrolling on button 4 and 5                       */
/*            scrolling width depend on shift and control keys     */
/*            without pressed key the scrollwidth 1/2 page         */
/*            with the control key                1   page         */
/*            with the shift key                  1   line         */
/*            Control + Shift is handled as Shift                  */
/*                                                                 */
/* INPUT:     Widget w              not relevant                   */
/*            XtPointer client_data really the scrollbar widget    */
/*            Xevent                the mosuse button event        */
/*                                                                 */
/* OUTPUT:    -                                                    */ 
/*                                                                 */ 
/* RETURN:    -                                                    */ 
/*                                                                 */ 
/* REMARKS:   Scrolling don't modify the selection                 */ 
/*                                                                 */ 
/*******************************************************************/ 

static void mouseScroll(Widget, XtPointer client_data, XEvent *event)
{
   Widget sb = (Widget)client_data;
   int value_return			  = 0;
   int slider_size_return	  = 0;
   int increment_return 	  = 0;
   int page_increment_return = 0;
	int count;
   int step;

   /* get a few value regarding the scrollbar conf. */
   XmScrollBarGetValues (sb, &value_return, &slider_size_return,
                        &increment_return, &page_increment_return);

   /* calculate the step wide according to the pressed keys */
   if ( event->xbutton.state & ShiftMask )
	{
	   step = 1;
	}
	else if ( event->xbutton.state & ControlMask )
	{
	   step = page_increment_return;
	}
	else
	{
	   step = page_increment_return >> 1;
	}
   
   if ( event->xbutton.button == Button4 )
	{
		value_return -= step;
	   if ( value_return < 0 )
		   value_return = 0;
	}
	else if ( event->xbutton.button == Button5 )
	{
      /* and the max value for increment */
      XtVaGetValues(sb, XmNmaximum, &count, NULL);
		value_return += step;
		if ( value_return > count - slider_size_return )
		   value_return = count - slider_size_return;
	}

   /* finally perform scrolling with the calculated step */
   if ( event->xbutton.button == Button4 ||
	     event->xbutton.button == Button5    )
	{
      XmScrollBarSetValues (sb, value_return, slider_size_return,
                        increment_return, page_increment_return, True);
   }
}

/*******************************************************************/ 
/*                                                                 */ 
/* NAME:      xmAddMouseEventHandler                               */ 
/*                                                                 */ 
/* FUNCTION:  Register the event handler for the button 4 and 5    */
/*                                                                 */
/* INPUT:     Widget w  The list/text widget                       */
/*                                                                 */
/* OUTPUT:    -                                                    */ 
/*                                                                 */ 
/* RETURN:    -                                                    */ 
/*                                                                 */ 
/* REMARKS:   -                                                    */ 
/*                                                                 */ 
/*******************************************************************/ 

void xmAddMouseEventHandler(Widget w)
{
   Widget wid;

   /* we need to pass the scrollbar widget to the handler */	
   XtVaGetValues(XtParent(w),XmNverticalScrollBar, &wid, NULL);

   /* handler for the scrolledList/ScrolledText */
   XtAddEventHandler(w, ButtonReleaseMask, False,
                     (XtEventHandler) mouseScroll, wid);
	/* and for the scrollbar itself */    
	XtAddEventHandler(wid, ButtonReleaseMask, False,
						  (XtEventHandler) mouseScroll, wid);
} 

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
