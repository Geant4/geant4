#include "main.h"
#include "torus.h"
#include <signal.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

double cos(double x);
double sin(double x);

#define EPSILON 1e-8
#define NBCERCLES 36
#define NBPOINTS  36

#define RAYON 100

int NbAff = 0;
double *px,*py,*pz;

void alarme (int sig)
{
		printf("NbAff = %d\n",NbAff);
		NbAff = 0;
		alarm (10);
}


int main (int argc,char *argv[]) {

  if (argc > 1) {
    printf("%s\n",argv[0]);
    printf("%s\n",argv[1]);
    if (strcmp(argv[1],"BVM") == 0) { 
      printf("BVM only\n");
      BVM_ONLY = 1;
    }
  }	

  /* Init X11 for 24 bpp */

  init_x ();

  /* Set up the alarm */
  
  signal (SIGALRM,alarme);
  alarm (10);

  while (event_x()) {						
    refresh ();
    aff ();
  }

  for(;;);		
  /* Free */

  free (px);
  free (py);
  free (pz);

  /* Close X */
  
  close_x();
  
  /* End */

  exit (0);
}

void aff () {
		NbAff ++;
		XPutImage(dis,win,gc,xim,0,0,0,0,W,H);

}

void refresh () {
  int i,j;
  long int V,R ;
  double d;
  static int t=1;

  double zoom = 0.20 + 0.01*t;
  double theta = 0.15 + 0.05*t;

  Intersect Inter ;

  memset ( buffer , 0 , 4*W*H );
  
  for (j=0;j<H;j++)
    for (i=0;i<W;i++)
      {
	//fprintf(stderr,"%d %d\n",i,j);
	
	R = 0;

	Inter.dx = -1.0*cos(theta);
	Inter.dy = 0.0;
	Inter.dz = -1.0*sin(theta);
		
	Inter.x = (65.0)*cos(theta) + 0.20*(i - W/2)*sin(theta); //45.0
	Inter.y = 0.20*(j - H/2) ;
	Inter.z = (65.0)*sin(theta) + 0.20*(i - W/2)*(-cos(theta));
						  
	Inter.R0 = 20.0; //20.0 25
	Inter.R1 = 14.0;  //4.0 19.0
	Inter.phi = 0.0;
	Inter.deltaphi = M_PI;
	d = DistanceToTorus (&Inter);
		
	if (d < 0)
	  {	
	    V = 0;
	    R = 255 ;
	  } else {
	    V = d*4;//4 16 32
	  }
	if (V > 255) {
	  V = 255 ;
	}
	
	if (V < 0) V = 0 ;
	buffer[4*(i*H+j)+0] = V;
	buffer[4*(i*H+j)+1] = V;
	buffer[4*(i*H+j)+2] = V+R;
      }
  t ++;
}

