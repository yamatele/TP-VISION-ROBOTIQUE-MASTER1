
/* Transformation de Hough sans sharing */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define WIDTH  1024                     /* TAILLE MAXI - reglee automatiquement par la lecture des fichiers*/
#define HEIGHT 1024                     /* TAILLE MAXI - reglee automatiquement*/
#define RESERV WIDTH*HEIGHT
#define PI    3.1415926535
#define RADIUS 800


#define MAXPOPULATION 5000

typedef struct {int theta;int ro;int life;float fitness;} mouche;

int i,j,k,ic,jc;
int width, height,radius;
float a,b,c,d;

unsigned char entree[15+RESERV],sortie[RESERV],hough[628*RADIUS],grad[RESERV],contour[RESERV];
int ientree[WIDTH][HEIGHT];
float isortie[WIDTH][HEIGHT];

/* fonctions */

void input(void);
void initialise(void);
void process(void);
void output(void);
float randfloat0(void);
float randfloat1(void);
float fitness(int);

/* déclarations spécifiques à l'application */


int theta,ro;
float SIN[628],COS[628];
float igrad[WIDTH][HEIGHT]; /* image de norme du gradient */
int THRESHOLD;
int POPULATION,GENERATIONS,PMUT,PCROSS,SAMPLES;
float MUTRO,MUTTHETA;
mouche indiv[MAXPOPULATION];
int generation;
int nb_evaluations=0, nb_gradients=0, nb_solutions=0;
float lambda;
float maxfitness;
float exigence;

/**********************************************************/
void input()
{
int truc;
FILE *fp;

/* lire le fichier parametres sous NFS */

  if ((fp = fopen("setup.txt","rb")))
    {fscanf(fp,"seuil %d\n",&THRESHOLD);
     fscanf(fp,"population %d\n",&POPULATION);
     fscanf(fp,"generations %d\n",&GENERATIONS);
     fscanf(fp,"pmutation %d\n",&PMUT);
     fscanf(fp,"pcrossover %d\n",&PCROSS);
     fscanf(fp,"mutabilite ro %f\n",&MUTRO);
     fscanf(fp,"mutabilite theta %f\n",&MUTTHETA);
     fscanf(fp,"samples %d\n",&SAMPLES);
     fscanf(fp,"exigence %f\n",&exigence);
     fclose(fp);

     printf ("seuil = %d\n",THRESHOLD);
     printf ("population = %d\n",POPULATION);
     if (POPULATION > MAXPOPULATION)
       {printf("\n debordement, population trop grande");
        POPULATION =MAXPOPULATION;
        printf ("population = %d\n",POPULATION);
       }
     printf ("nombre de generations = %d\n",GENERATIONS);
     if (PMUT<0) PMUT=0;
     else if(PMUT>100) PMUT=100;
     if (PCROSS<0)PCROSS=0;
     else if(PCROSS>100)PCROSS=100;
     printf("prob.mut.=%d\tprob.cross.=%d\n",PMUT,PCROSS);
     printf("mutabilites: ro=%f\t theta=%f\n",MUTRO,MUTTHETA);
     printf("nombre d'echantillons = %d\n",SAMPLES);
     printf("coeff d'exigence = %f\n",exigence);

    }
  else {printf("fichier setup non trouve\n"); exit(1);}

/* lire le fichier image sous NFS */

  if ((fp = fopen("entree.pbm","rb")))
    {fscanf(fp,"P%d\n",&truc); /* truc = 4, 5 ou 6 */
     fscanf(fp,"%d %d\n",&width,&height);
     fclose (fp);
     if(width*height > WIDTH*HEIGHT) {printf("format image trop grand\n");exit(1);}
     printf("width = %d  height = %d \n", width, height);
    }
  else {printf("fichier entree.pbm non trouve\n"); exit(1);}
  
  radius = height*height+width*width;
  radius = sqrt(radius)/2;


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
      grad[i+width*j] = (unsigned char) (igrad[i][j]); 

for (i=0;i<width;i++)  for (j=0;j<height;j++)
      sortie[i+width*j] = (unsigned char) (isortie[i][j]); 


/* ecrire les fichiers sous NFS */

fp=fopen("grad.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(grad,1,width*height,fp);
fclose(fp);


fp=fopen("sortie.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(sortie,1,width*height,fp);
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

/* initialise la stratégie d'évolution */
for(i=0;i<POPULATION;i++)
  {
   indiv[i].theta=628*randfloat0();
   indiv[i].ro=radius*randfloat0();
   indiv[i].fitness=0;
   indiv[i].life=1;
                      /* life = 0 vivant évalué (on a calculŽ sa fitness)
                                1 vivant non évalué (on n'a pas encore calculŽ sa fitness)
                                2 mort à remplacer par mutation
                                3 mort à remplacer par croisement
                      */
  }
maxfitness=0.;
}

/*******************************************/

void process()

{

 for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
   {
    a = abs(ientree[i][j+1] - ientree[i][j-1]);
    c = abs(ientree[i+1][j] - ientree[i-1][j]);
    igrad[i][j] = 2*(a+c);
   }

printf("\ngradient calcule");

/* maintenant on change de methode */

maxfitness=0;
for(generation=0;generation<GENERATIONS;generation++)
  {


           /* mise à jour de la fitness */
//life = 0 individu vivant ŽvaluŽ
//life = 1 individu vivant non ŽvaluŽ

   for(i=0;i<POPULATION;i++)
     {
      if(indiv[i].life == 1) {indiv[i].fitness=fitness(i);indiv[i].life=0;}
     }
   printf("\nmaxfitness %5.0f",maxfitness);
           /* tournois de sélection */
   for(k=0;k<(POPULATION*PMUT)/100;k++)
     {
      i=POPULATION*randfloat0();
      j=POPULATION*randfloat0();
      if(indiv[i].fitness<indiv[j].fitness) indiv[i].life=2;else indiv[j].life=2; /* je tue le plus mauvais */
     }
   for(k=0;k<(POPULATION*PCROSS)/100;k++)
     {
      i=POPULATION*randfloat0();
      j=POPULATION*randfloat0();
      if(indiv[i].fitness<indiv[j].fitness) indiv[i].life=3;else indiv[j].life=3; /* je tue le plus mauvais */

     }
           /* activation des opérateurs génétiques */
    for(i=0;i<POPULATION;i++)
     {
      if(indiv[i].life==2)   /* individu a remplacer par mutation */
       {
        j=POPULATION*randfloat0();                          /* choix du parent j */
        a=fabsf(indiv[j].ro+MUTRO*randfloat1());
        if(a<1)a=1;else if(a>radius)a=radius-1;
        indiv[i].ro=a;

        a=indiv[j].theta+MUTTHETA*randfloat1();
        if(a<0) a+=628; else if (a>628) a -= 628;
        indiv[i].theta=a;

        indiv[i].life=1;
       }
      else if (indiv[i].life==3) /*individu a remplacer par croisement */
       {
        j=POPULATION*randfloat0();
        k=POPULATION*randfloat0();                      /* choix des 2 parents j et k */
        lambda=-.5+2*randfloat0();

        a=fabsf(lambda*indiv[j].ro+(1-lambda)*indiv[k].ro);
        if(a<1)a=1;else if(a>radius)a=radius-1;
        indiv[i].ro=a;

        a=lambda*indiv[j].theta+(1-lambda)*indiv[k].theta;
        if(a<0) a+=628; else if (a>628) a -= 628;
        indiv[i].theta=a;

        indiv[i].life=1;
       }
     }
  }


/*maintenant on a une population, quoi en faire?  */

for (i=1;i<width-1;i++) for (j=1;j<height-1;j++) isortie[i][j]=ientree[i][j]/2;

for(k=0;k<POPULATION;k++)
 if(indiv[k].fitness > exigence * maxfitness )
  {
   nb_solutions ++;
   printf("\n une solution, theta=%d,ro=%d,qualite=%f",indiv[k].theta,indiv[k].ro,indiv[k].fitness/maxfitness);
   a=COS[indiv[k].theta];b=SIN[indiv[k].theta];
   for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
    {
     ic=i-width/2;jc=j-height/2;    /* coordonnees centrées */
     if (fabsf(indiv[k].ro-ic*a-jc*b)<.5) isortie[i][j]=255;

    }
  }
 if(nb_solutions<2)
for(k=0;k<POPULATION;k++)
 if(indiv[k].fitness > exigence * exigence*maxfitness )
  {
   printf("\n une solution, theta=%d,ro=%d,qualite=%f",indiv[k].theta,indiv[k].ro,indiv[k].fitness/maxfitness);
   a=COS[indiv[k].theta];b=SIN[indiv[k].theta];
   for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
    {
     ic=i-width/2;jc=j-height/2;    /* coordonnees centrées */
     if (fabsf(indiv[k].ro-ic*a-jc*b)<.5) isortie[i][j]=255;

    }
  }
   printf("\nTermine: %d solutions, %d evaluations, %d gradients calcules \n\n",nb_solutions,nb_evaluations,nb_gradients);
}
/*******************************************/
float fitness(int i)
{
 float mu;int x,y,xc,yc;
 c=0;

 theta=indiv[i].theta;
 a=COS[theta];b=SIN[theta];
 ro=indiv[i].ro;
/*        printf("\ntheta=%d,ro=%d",theta,ro); */

 for(j=0;j<SAMPLES;j++)
  {
   mu=radius*randfloat1();

   xc=ro*a+mu*b;
   yc=ro*b-mu*a;
   x=xc+width/2;
   y=yc+height/2;
   if((x>0)&&(y>0)&&(x<width)&&(y<height))
     {nb_gradients++;
      if(igrad[x][y]>THRESHOLD) c += igrad[x][y]; //on pourrait essayer c += 1;
     }
  }

if(c>maxfitness)maxfitness=c;
nb_evaluations ++;
return c;
}

/*******************************************/
float randfloat0(void)
{
  float alea, numerator, denominator;

  numerator = rand() & 0x7fff;
  denominator = 0x7fff;

  alea = numerator/denominator;
  return alea;  /* equireparti sur [0,1] */
}
/**************************************************************** RANDFLOAT1 */
float randfloat1(void)
{
  float alea, numerator, denominator;

  numerator = rand() & 0x7fff;
  denominator = 0x7fff;

  alea = 2. *(numerator/denominator) - 1.;
  alea = alea * (.2 + .8 * alea * alea ) ;
  return alea;  /* reparti sur [-1,1], concentre pres de 0 */
}

/*******************************************/

int main()
{
input();
initialise();
process();
output();
}
