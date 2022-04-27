
/* Transformation de Hough version 2016 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define WIDTH  2048                     /* TAILLE MAXI - reglee automatiquement par la lecture des fichiers*/
#define HEIGHT 2048                     /* TAILLE MAXI - reglee automatiquement*/
#define RESERV WIDTH*HEIGHT
#define PI    3.1415926535
#define RADIUS 1500

int i,j,ic,jc;
int width, height;
int seuil;
float a,b,c,d;
int THRESHOLD;

unsigned char entree[15+RESERV],sortie[RESERV],hough[628*RADIUS],grad[RESERV],contour[RESERV];
int ientree[WIDTH][HEIGHT];
float isortie[WIDTH][HEIGHT];

/* fonctions */

void input(void);
void initialise(void);
void process(void);
void output(void);

/* déclarations spécifiques à l'application */


int theta,ro;
float SIN[628],COS[628];
float ihough[628][RADIUS];
float igrad[WIDTH][HEIGHT]; /* image de norme du gradient */
float icontour[WIDTH][HEIGHT]; /* image de contours */


/**********************************************************/
void input()
{
int truc;
FILE *fp;

/* lire le fichier parametres sous NFS */

/*
images .pbm monochrome = .pgm
P5
800 700
255
piwehoejgojg;skj;werjqoiwlshbsldkjn
*/

  if ((fp = fopen("setup.txt","rb")))
    {fscanf(fp,"seuil %d\n",&THRESHOLD);
     fclose(fp);
     printf ("seuil = %d\n",THRESHOLD);
    }
  else {printf("fichier setup non trouve'\n"); exit(1);}

/* lire le fichier image sous NFS */

  if ((fp = fopen("entree.pbm","rb")))
    {fscanf(fp,"P%d\n",&truc); /* truc = 4, 5 ou 6 */
     fscanf(fp,"%d %d\n",&width,&height);
     fclose (fp);
     if(width*height > WIDTH*HEIGHT) {printf("format image trop grand\n");exit(1);}
     printf("width = %d  height = %d \n", width, height);
    }
  else {printf("fichier entree.pbm non trouve\n"); exit(1);}
  

  fp = fopen("entree.pbm","rb");
  fread (entree,1,RESERV+15,fp);
  fclose(fp);

  if (truc != 5)
    {printf("\n attention ce programme a besoin d'images pbm monochrome P5 en entree - abandon\n");
     exit(1);
    }

/* lecture des parametres et codage en entier */

 for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      ientree[i][j] = (int) entree[15+(i+width*j)];
  
}
/*******************************************/
void output()
{
FILE *fp;

/* transformer les tableaux en fichiers .pbm monochromes */

for (i=0;i<width;i++)  for (j=0;j<height;j++)
      sortie[i+width*j] = (unsigned char) (isortie[i][j]); 


for (i=0;i<width;i++)  for (j=0;j<height;j++)
      grad[i+width*j] = (unsigned char) (igrad[i][j]); 


for (i=0;i<width;i++)  for (j=0;j<height;j++)
      contour[i+width*j] = (unsigned char) (icontour[i][j]); 


for(theta=0;theta<628;theta++)  for(ro=0;ro<RADIUS;ro++)
      hough[theta+628*ro] = (unsigned char) (ihough[theta][ro]);

/* ecrire les fichiers sous NFS */

fp=fopen("sortie.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(sortie,1,width*height,fp);
fclose(fp);

fp=fopen("grad.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(grad,1,width*height,fp);
fclose(fp);


fp=fopen("contour.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(contour,1,width*height,fp);
fclose(fp);

fp=fopen("hough.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",628,RADIUS);
fwrite(hough,1,628*RADIUS,fp);
fclose(fp);

}
/*******************************************/
void initialise()   
{
/* tabule les fonctions trigo */
for(theta=0;theta<628;theta++)
 {
  COS[theta]=cos(theta/100.);
  SIN[theta]=sin(theta/100.);
 }

/*nettoie l'accumulateur*/
for(theta=0;theta<628;theta++)
  {for(ro=0;ro<RADIUS;ro++)
    ihough[theta][ro]=0.;
  }
}

/*******************************************/


/*******************************************/

void process()

{
int orientation;   /* orientation du gradient */
float maxhough = 0.;
float truc;
int thetamax,romax;

 for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
   {
    a = abs(ientree[i][j+1] - ientree[i][j-1]);
    b = abs(ientree[i+1][j+1]-ientree[i-1][j-1]);
    c = abs(ientree[i+1][j] - ientree[i-1][j]);
    d = abs(ientree[i+1][j-1]-ientree[i-1][j+1]); 
    igrad[i][j] = a+b+c+d;

    if(a>b)
      {if(a>c)
         {if(a>d)  orientation=1; /* a est dominant */
          else     orientation=4; /* d est dominant */
         }
       else {if (c<d) orientation=4; /* d dominant */
             else orientation=3; /* c dominant */
            }
      }
     else  /*a est hors jeu */
       {if (b>c)
          {if(b>d) orientation=2; /* b dominant */
           else orientation = 4; /* d dominant */
          }
        else {if (c<d) orientation=4; /* d dominant */
              else orientation=3; /* c dominant */
             }
       }


/* on connait l'orientation du gradient: cherchons maintenant les points de contour */

   if(igrad[i][j]>THRESHOLD)
      {if (orientation==1) {if (icontour[i][j-1]==0) icontour[i][j]=255;}
       else if (orientation==2) {if (icontour[i-1][j-1]==0) icontour[i][j]=255;}
       else if (orientation==3) {if (icontour[i-1][j]==0) icontour[i][j]=255;}
       else                   {if (icontour[i-1][j+1]==0) icontour[i][j]=255;}
      }
   }
printf("\ncontours calcules");

/* maintenant remplissons l'accumulateur de Hough à partir des points de contour détectés */

 for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
   {
   
   
   
   for (theta=0;theta<628;theta++)
    {
     ro = (int)(ic*COS[theta]+jc*SIN[theta]);
     if(ro>0) if(ro<RADIUS-2)
        {
        truc=ihough[theta][ro];
        ihough[theta][ro] += igrad[i][j]/1000;
        if(maxhough<=truc) {maxhough=truc+1;thetamax=theta;romax=ro;}
          
        }
    }
    
    /*
    if(icontour[i][j]==255)
      {
       ic=i-width/2;jc=j-height/2;    // coordonnees centrees
       for(theta=0;theta<628;theta++)
         {
          ro = (int)(ic*COS[theta]+jc*SIN[theta]);
          if(ro>0)  // attention à virer les valeurs négatives de ro, une sur deux en gros...
           if(ro<RADIUS-2)
            {
             truc=ihough[theta][ro];
             ihough[theta][ro] = truc+.5;
             ihough[theta][ro+1] += .25;
             ihough[theta][ro-1] += .25;

            }
          if(maxhough<=truc) {maxhough=truc+1;thetamax=theta;romax=ro;}
         }
      }
      */
   }
 printf("\nHough: max = %f\tthetamax=%d\tromax=%d\n",maxhough,thetamax,romax);
ihough[thetamax][romax]=255;

/* image de sortie = image d'entrée avec la droite détectée en surimpression */


 a=COS[thetamax];b=SIN[thetamax];
 for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
  {
   ic=i-width/2;jc=j-height/2;    /* coordonnees centrées */
   if (fabsf(romax-ic*a-jc*b)<1) isortie[i][j]=255;
   else isortie[i][j]=ientree[i][j];
  }
 printf("Hough: image de sortie calculee, travail termine\n\n");
}

/*******************************************/

int main()
{
input();
initialise();
    printf("this is the classical Hough transform\n");
process();
output();
}
