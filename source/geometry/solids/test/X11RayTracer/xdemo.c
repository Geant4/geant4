#include "xdemo.h"

void init_x () {

  dis = XOpenDisplay  ( NULL );
  if ( dis == NULL)
    {
      fprintf ( stderr , "Display: connection impossible\n");
      exit (1);
    }

  screen = DefaultScreen ( dis );
  depth  = DefaultDepth  ( dis , screen );
  width  = DisplayWidth  ( dis , screen );
  height = DisplayHeight ( dis , screen );

  if ( depth != 24 )
    {
      fprintf ( stderr , "This code works only in 24 bpp..\n" );
      XCloseDisplay  ( dis );
      exit (1);
    }

  winRoot                  = DefaultRootWindow ( dis );
  winAttr.border_pixel     = BlackPixel ( dis , screen );
  winAttr.background_pixel = BlackPixel ( dis , screen );
  winMask                  = CWBackPixel | CWBorderPixel;
 
  win = XCreateWindow ( dis , winRoot , X , Y , W , H , BORDER , depth ,
			InputOutput , CopyFromParent , winMask , &winAttr );

  XStoreName ( dis , win , NAME );

  XSelectInput ( dis , win , KeyPressMask );

  winHint.flags                           = PPosition | PMinSize | PMaxSize ;
  winHint.x                               = X;
  winHint.y                               = Y;
  winHint.max_width  = winHint.min_width  = W;
  winHint.max_height = winHint.min_height = H;
  XSetWMNormalHints ( dis , win , &winHint );

  XClearWindow ( dis , win );
  XMapRaised ( dis , win );
  XFlush ( dis );

  buffer = ( char * ) malloc ( 4*W*H );
  xim = XCreateImage ( dis , CopyFromParent , depth , ZPixmap , 0 ,
		       buffer , W , H , 32 , W*4 );

  gcVal.foreground = 0;
  gcVal.background = 0;
  gcMask           = GCForeground | GCBackground;

  gc = XCreateGC ( dis , win , gcMask , &gcVal );

}


int event_x () {

  XEvent XEv;

  return ( ! XCheckWindowEvent ( dis , win , KeyPressMask , &XEv ) ); 

}


void close_x () {

  XFreeGC        ( dis , gc );
  XDestroyImage  ( xim );
  XDestroyWindow ( dis , win );
  XCloseDisplay  ( dis );     

};

