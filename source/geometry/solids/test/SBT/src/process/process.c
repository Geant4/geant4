/*
  MEDERNACH Emmanuel
  Aug  8, 2000
 */

 #include <stdio.h>
 #include <unistd.h>
 #include <sys/types.h>
 #include <sys/wait.h>

 /*
 TODO:
 Add support for 'error management' 
  */

 pid_t wait(int *status) ;

 void read_choice(int *choice)
 {

   if (choice == NULL) return ;
   (*choice) = -1;

   printf("Possible volumes are :\n");
   printf("1: BOX\n");
   printf("2: TUBS\n");
   printf("3: CONE\n");
   printf("4: PCON\n");
   printf("5: PGON\n");
   printf("6: PGON2\n");
   printf("7: CONE2 (Emm)\n");
   printf("8: THINBOX (Emm)\n");
   while (((*choice) < 1) || ((*choice) > 8))
     {
       scanf("%d",choice);
     }
 }

 void write_ini(const char *init,int choice,double x,double y,double z,double dx,double dy,double dz)
 {
   FILE *fi= NULL;

   fi = fopen(init,"w");
   if (fi == NULL)
     {
       perror("Some problems to write in the Init.fred file");
       exit (1);
     }

   switch (choice)
     {
     case 1:
       fprintf(fi,"/fred/volume BOX\n");
       break;
     case 2:
       fprintf(fi,"/fred/volume TUBS\n");
       break;
     case 3:
       fprintf(fi,"/fred/volume CONE\n");
       break;
     case 4:
       fprintf(fi,"/fred/volume PCON\n");
       break;
     case 5:
       fprintf(fi,"/fred/volume PGON\n");
       break;
     case 6:
       fprintf(fi,"/fred/volume PGON2\n");
       break;
     case 7:
       fprintf(fi,"/fred/volume CONE2\n");
       break;
     case 8:
       fprintf(fi,"/fred/volume THINBOX\n");
       break;
     default:
       fprintf(stderr,"%d is not between 1 and 8\n",choice);
       exit(2);
       break;
     }
   fprintf(fi,"/fred/gun G4\n");
   fprintf(fi,"/gun/position %f %f %f\n",x,y,z);
   fprintf(fi,"/gun/direction %f %f %f\n",dx,dy,dz);
   fprintf(fi,"/run/initialize \n");
   fprintf(fi,"/run/beamOn\n");
   fclose(fi);
 }

 int main ()
 {
   pid_t son = -1;
   int status = 0;
   char init[]="Init.fred";
   char *param [3];

   int choice = -1 ;
   double posx,posy,posz;
   double dirx,diry,dirz;

   double startx,starty,startz;
   double beginx,beginy,beginz;
   double endx,endy,endz;
   double deltax,deltay,deltaz;
   double begindirx,enddirx,stepdirx;
   double begindiry,enddiry,stepdiry;
   double begindirz,enddirz,stepdirz;

   int TakeOtherPoint = 0;


   /* write an init file */

   read_choice(&choice);


   posx = posy = posz = 0;
   dirx = diry = dirz = 1;


   printf("You could choose a delta <= 0 to have only a plane\n");
   /* rotation in a cone for example could be performed */
   /* to have only a y=0 plane to search */ 
   printf("Choose a beginning, an end and a step for x (0 200 10)\n");
   scanf("%lf %lf %lf",&beginx,&endx,&deltax);
   printf("Choose a beginning, an end and a step for y (0 200 10)\n");
   scanf("%lf %lf %lf",&beginy,&endy,&deltay);
   printf("Choose a beginning, an end and a step for z (0 200 10)\n");
   scanf("%lf %lf %lf",&beginz,&endz,&deltaz);

   /* you could want to continue after a previous crash 
      (X could hangs if swap is full)
   */

   printf("you could want to continue after a previous crash\n");
   printf("In this case just indicate from where do you want to start (x y z)\n\n");
   scanf("%lf %lf %lf",&startx,&starty,&startz);

   printf("you could now have a grid of direction to point\n");
   printf("Choose a beginning, an end and a step for all direction (-30 30 1)\n");
   printf("\nfor direction x : ");
   scanf("%lf %lf %lf",&begindirx,&enddirx,&stepdirx);
   printf("\nfor direction y : ");
   scanf("%lf %lf %lf",&begindiry,&enddiry,&stepdiry);
   printf("\nfor direction z : ");
   scanf("%lf %lf %lf",&begindirz,&enddirz,&stepdirz);

   /* we want also != step for each direction */

   posx = startx ;
   posy = starty ;
   posz = startz ;


   printf("We start at %f %f %f\n",posx,posy,posz);

   while (posx < endx)
     { /* loop on x */

       printf("Posx = %f\n",posx);
       while (posy < endy)
	 { /* loop on y */

	   printf("Posy = %f\n",posy);
	   while (posz < endz)
	     { /* loop on z */

	       printf("\tPosz = %f\n",posz);
	       TakeOtherPoint = 0;
	       for (dirx=begindirx;(dirx<enddirx) && (TakeOtherPoint == 0);dirx+=stepdirx)
		 for (diry=begindiry;(diry<enddiry) && (TakeOtherPoint == 0);diry+=stepdiry)
		   for (dirz=begindirz;(dirz<enddirz) && (TakeOtherPoint == 0);dirz+=stepdirz)
		     {
		       write_ini(init,choice,posx,posy,posz,dirx,diry,dirz);

		       /* call to Fred */
		       son = fork();
		       switch(son)
			 {
			 case -1:
			   perror("Fork process error");
			   exit(1);
			   break;
			 case 0:
			   /* execvp (Emm)fred with init file */
			   param[0] = "fred";
			   param[1] = init;
			   param[2] = NULL;

			   /* no output nor input */
			   close (0);
			   close (1);
			   close (2);

			   if  (execvp("fred",param) == -1)
			     {
			       fprintf(stderr,"fred program not found ..\n");
			       exit(1);
			     }
			   break;
			 default:
			   /* wait for the other process */
			   waitpid(son,&status,0);
			   if (WIFEXITED(status) == 0)
			     {
			       /* some errors */
			       /* keep the current init */
			       fprintf(stderr,"I catch an error : bad return value %f %f %f %f %f %f\n",posx,posy,posz,dirx,diry,dirz);
			       /* we take another point */
			       /* Ok this is not really the way to do ..*/
			       TakeOtherPoint = 1;
			     }
			   else
			     if (WIFSIGNALED(status))
			       {
				 /* some errors */
				 /* keep the current init */
				 fprintf(stderr,"I catch an error : Signal terminated %f %f %f %f %f %f\n",posx,posy,posz,dirx,diry,dirz);
				 /* we take another point */
				 /* Ok this is not really the way to do ..*/
				 TakeOtherPoint = 1;
			       }
			   break;
			 }
		     }
	       /* loop on z */
	       posz += deltaz ;
	     }
	   posz = beginz ;

	   /* loop on y */
	   posy += deltay ;
	 }
       posy = beginy ;

       /* loop on x */
       posx += deltax ;
     }
   posx = beginx ;

 }

