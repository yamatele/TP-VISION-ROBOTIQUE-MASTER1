/* Evolutionary circle detection */
// bangorcircle01.c inspired from evocircle26.c
// with setup file renamed circlesetup.txt
// bangorcircle 02 introduces minradius into the new circlesetup.txt file
// and a more specific fitness function
//bangorcircle 03 doesnt use gradient anymore as the objects in original image already have bright contoursot
//bangorcircle05 the fitness parameters are now set in the setup file. Added a "darkness" factor.
// newcircle05.c pour les TP EFREI
//06.c ah=jout d'un critere de gradient sur le bord

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define WIDTH  3000                    /* TAILLE MAXI - reglee automatiquement par la lecture des fichiers*/
#define HEIGHT 3000                     /* TAILLE MAXI - reglee automatiquement*/
#define RESERV WIDTH*HEIGHT
#define MAXPOPULATION 1000000

////////////////////////////////
typedef struct {int xc; int yc; int r; float fitness; int retained;} object;
///////////////////////////////
int x,y,xc,yc,r,theta, i,j,k,l,truc;
int width, height;
float a,b,c,d;
float meangreylevel=0.;   //not used in this version!
int SHARINGRADIUS;

unsigned char input[15+RESERV],output[RESERV],hough[RESERV*20],grad[RESERV],contour[RESERV];
int iinput[WIDTH][HEIGHT];
float ioutput[WIDTH][HEIGHT];
float tabsharing[2+WIDTH][2+HEIGHT];

/* functions */

void readinput(void);
void initialise(void);
void process(void);
void showcircles(void);
void writeoutput(void);
float randfloat0(void);
float randfloat1(void);
float fitness(int);

/* application-specific declarations */

float SIN[628],COS[628];
float igrad[WIDTH][HEIGHT]; /* gradient norm */
int POPULATION,GENERATIONS,PMUTPOS,PMUTRAD,PCROSS,PIMM,SAMPLES,MINRADIUS,MAXRADIUS;
float edge, contrast, homogeneity, darkness, prefdiameter, border;

float MUTPOSITION, MUTRADIUS;
object indiv[MAXPOPULATION];
int generation;
int nb_evaluations=0, nb_solutions=0;
float lambda;
float maxfitness,maxshfitness;
float exigence;int numberofcircles;

/**********************************************************/
void readinput()
{
int truc;
FILE *fp;

/* read setup parameter file under NFS */

  if ((fp = fopen("circlesetup.txt","rb")))
    {
     fscanf(fp,"population %d\n",&POPULATION);  //10000
     fscanf(fp,"generations %d\n",&GENERATIONS);  //1000000
	 fscanf(fp,"pimmigration %d\n",&PIMM);  //1
     fscanf(fp,"pmutationposition %d\n",&PMUTPOS);//3
     fscanf(fp,"pmutationradius %d\n",&PMUTRAD);//3
     fscanf(fp,"pcrossover %d\n",&PCROSS);  //1
     fscanf(fp,"mutability of position %f\n",&MUTPOSITION); // 20
     fscanf(fp,"mutability of radius %f\n",&MUTRADIUS);  //  10
     fscanf(fp,"minradius %d\n",&MINRADIUS); //5
     fscanf(fp,"maxradius %d\n",&MAXRADIUS); //50
     fscanf(fp,"samples %d\n",&SAMPLES); //30
     fscanf(fp,"number of circles %d\n", &numberofcircles); //75
     fscanf(fp,"sharing radius %d\n", &SHARINGRADIUS); //25
     fscanf(fp,"edge %f\n",&edge); //3
     fscanf(fp,"contrast %f\n",&contrast); //24
     fscanf(fp,"homogeneity %f\n",&homogeneity); //8
     fscanf(fp,"darkness %f\n", &darkness); //8
     fscanf(fp,"prefdiameter %f\n",&prefdiameter); //4
     fscanf(fp,"border contrast %f\n", &border); //4
     fclose(fp);

	
     printf ("population = %d\n",POPULATION);
     if (POPULATION > MAXPOPULATION)
       {printf("\n population overflow");
        POPULATION =MAXPOPULATION;
        printf ("population = %d\n",POPULATION);
       }
     printf ("number of generations = %d\n",GENERATIONS);
     if (PMUTPOS<0) PMUTPOS=0;
     else if(PMUTPOS>100) PMUTPOS=100;
		if (PMUTRAD<0) PMUTRAD=0;
		else if(PMUTRAD>100) PMUTRAD=100;
		
		
     if (PCROSS<0)PCROSS=0;
     else if(PCROSS>100)PCROSS=100;
     if (PIMM<0) PIMM=0;
     else if(PIMM>100) PIMM=100;
     printf("Prob. imm.=%d\tprob.mut.pos.=%d\tprob.mut.rad.=%d\tprob.cross.=%d\n",PIMM,PMUTPOS,PMUTRAD,PCROSS);
     printf("mutabilities: position=%f\t radius=%f\n",MUTPOSITION,MUTRADIUS);
     printf("minimum circle radius = %d\n", MINRADIUS);
     printf("maximum circle radius = %d\n", MAXRADIUS);
     printf("number of samples = %d\n",SAMPLES);
     printf("number of circles = %d\n",numberofcircles);    //now useless?
     printf("sharingradius = %d\n", SHARINGRADIUS);

    }
  else {printf("file circlesetup not found\n"); exit(1);}

/* read image file under NFS */

  if ((fp = fopen("input.pbm","rb")))
    {fscanf(fp,"P%d\n",&truc); /* truc = 4, 5 ou 6 */
     fscanf(fp,"%d %d\n",&width,&height);
     fclose (fp);
     if(width*height > WIDTH*HEIGHT) {printf("image too large\n");exit(1);}
     printf("width = %d  height = %d \n", width, height);
    }
  else {printf("file input.pbm not found\n"); exit(1);}
  

  fp = fopen("input.pbm","rb");
  fread (input,1,RESERV+15,fp);
  fclose(fp);

  if (truc != 5)
    {printf("\n this program needs pbm P5 input images - abandon\n");
     exit(1);
    }

/* reading parameters and coding them as integers */

 for (i=0;i<width;i++) for (j=0;j<height;j++)
      iinput[i][j] = (int) input[15+(i+width*j)];
  
}
/*******************************************/
void writeoutput()
{
FILE *fp;

/* transform arrays into grey level .pbm P5 files */

for (i=0;i<width;i++)  for (j=0;j<height;j++)
      grad[i+width*j] = (unsigned char) (igrad[i][j]/8); 

for (i=0;i<width;i++)  for (j=0;j<height;j++)
      output[i+width*j] = (unsigned char) (ioutput[i][j]); 

/* create files under NFS */
fp=fopen("grad.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(grad,1,width*height,fp);
fclose(fp);

fp=fopen("output.pbm","wb");
fprintf(fp,"P5\n%d %d\n255\n",width,height);
fwrite(output,1,width*height,fp);
fclose(fp);

}
/*******************************************/
void initialise()   
{
//calculate mean grey level of image
for(i=0;i<width;i++) for (j=0;j<height;j++) meangreylevel += iinput[i][j];
 meangreylevel /= width*height;                                //just a counter
x=indiv[i].xc; y=i;
printf("mean grey level = %f\n",meangreylevel);

//calculate gradient image

//for (i=0;i<width;i++) for (j=0;j<height;j++) igrad[i][j]=iinput[i][j]; //not calculating grad image!

for (i=0;i<width;i++) for (j=0;j<height;j++) igrad[i][j]=0;
	
for (i=1;i<width-1;i++) for (j=1;j<height-1;j++)
   {
    a=abs(iinput[i+1][j-1]+iinput[i+1][j]+iinput[i+1][j+1]-iinput[i-1][j-1]-iinput[i-1][j]-iinput[i-1][j+1]);
    b=abs(iinput[i+1][j-1]+iinput[i][j-1]+iinput[i-1][j-1]-iinput[i+1][j+1]-iinput[i][j+1]-iinput[i-1][j+1]);
    c=(a+b);
    igrad[i][j] += 2*c;
    igrad[i-1][j] += c;
    igrad[i+1][j] += c;
    igrad[i][j-1] += c;
    igrad[i][j+1] += c;
   }
printf("\ngradient calculated\n");

/* tabulate trigo functions */
for(i=0;i<628;i++)
 {COS[i]=cos(i/100.);SIN[i]=sin(i/100.);}
	
   maxfitness = -10000.;
	
	/*initialise population!*/
	
	for(i=0;i<POPULATION;i++)
	  {
       c=0;
	   while((c<MINRADIUS)||(c>MAXRADIUS)) c = MINRADIUS+(MAXRADIUS-MINRADIUS)*sqrt(randfloat0());
       indiv[i].r=c;
       indiv[i].xc=MAXRADIUS +(width-2.*MAXRADIUS)*randfloat0();
       indiv[i].yc=MAXRADIUS +(height-2.*MAXRADIUS)*randfloat0();
       tabsharing[(int)indiv[i].xc/SHARINGRADIUS][(int)indiv[i].yc/SHARINGRADIUS] += 1;
	   indiv[i].fitness=fitness(i);
	   indiv[i].retained = 0;
      }	
	printf("initialisation successful\n");
}

/*******************************************/
void process()
{
    float averagefitness, coeff, coefficient;
	int counter;
	 FILE *fp;
	 int radius[POPULATION], stat[100];
	 for(k=0;k<100;k++) stat[k]=0;
 ///////////////BEGIN GENERATIONS/////////////////	
for(generation=0;generation<=GENERATIONS;generation++)
  {
  //printf("\ngen=%d",generation);
//replace by immigration

	  for(truc=0;truc<PIMM;truc++)
	  {
		  i=(POPULATION-1)*randfloat0();
		  j=(POPULATION-1)*randfloat0();
		  if(j!=i)
		  {
			if(indiv[i].fitness /(.01+tabsharing[(int)(indiv[i].xc/SHARINGRADIUS)][(int)(indiv[i].yc/SHARINGRADIUS)])
			   <indiv[j].fitness /(.01+tabsharing[(int)(indiv[j].xc/SHARINGRADIUS)][(int)(indiv[j].yc/SHARINGRADIUS)]))
			k=i;else k=j;
	//we kill k and replace it with a brand new individual
			tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] -= 1;
			c=0; while ((c<MINRADIUS)||(c>MAXRADIUS))c=MINRADIUS+(MAXRADIUS-MINRADIUS)*sqrt(randfloat0());
			indiv[k].r=c;
			indiv[k].xc=MAXRADIUS +(width-2.*MAXRADIUS)*randfloat0();
			indiv[k].yc=MAXRADIUS +(height-2.*MAXRADIUS)*randfloat0();
			tabsharing[(int)indiv[k].xc/SHARINGRADIUS][(int)indiv[k].yc/SHARINGRADIUS] += 1;
			indiv[k].fitness=fitness(k);
			if (maxfitness<indiv[k].fitness) maxfitness=indiv[k].fitness;
		  }
	  }
	  
// replace by mutation on centre position

   for(truc=0;truc<PMUTPOS;truc++)
     {
      i=(POPULATION-1)*randfloat0();
      j=(POPULATION-1)*randfloat0();

      if(j!=i)
        {
		 if(indiv[i].fitness /(.01+tabsharing[(int)(indiv[i].xc/SHARINGRADIUS)][(int)(indiv[i].yc/SHARINGRADIUS)])
		  < indiv[j].fitness /(.01+tabsharing[(int)(indiv[j].xc/SHARINGRADIUS)][(int)(indiv[j].yc/SHARINGRADIUS)]))
		 {k=i;l=j;} else {k=j;l=i;}
                    //we call k the loser and l the winner
	     tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] -= 1;
          // choose parent l (winner) as gene giver and k as child

	     indiv[k].r=indiv[l].r;
	     c=indiv[l].r;
	     a=0; while ((a<c)||(a>width-c)) a=indiv[l].xc+MUTPOSITION*randfloat1();
         indiv[k].xc=a;
		 b=0; while ((b<c)||(b>height-c)) b=indiv[l].yc+MUTPOSITION*randfloat1();
		 indiv[k].yc=b;
         tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] += 1;
         indiv[k].fitness=fitness(k);
         if (maxfitness<indiv[k].fitness) maxfitness=indiv[k].fitness;
          }
     }
	  
	  // replace by mutation on radius

for(truc=0;truc<PMUTRAD;truc++)
	{
	i=(POPULATION-1)*randfloat0();
	j=(POPULATION-1)*randfloat0();

	if(j!=i)
		{
		if(indiv[i].fitness /(.01+tabsharing[(int)(indiv[i].xc/SHARINGRADIUS)][(int)(indiv[i].yc/SHARINGRADIUS)])
		   <indiv[j].fitness /(.01+tabsharing[(int)(indiv[j].xc/SHARINGRADIUS)][(int)(indiv[j].yc/SHARINGRADIUS)]))
		{k=i;l=j;} else {k=j;l=i;}
			  //we call k the loser and l the winner
		tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] -= 1;
			   // choose parent l (winner) as gene giver and k as child
		c=0;
		while ((c<MINRADIUS)||(c>MAXRADIUS)) c=(indiv[l].r+MUTRADIUS*randfloat1());
		indiv[k].r=c;
		indiv[k].xc=indiv[l].xc;
		indiv[k].yc=indiv[l].yc; //ne me sois point rebelle puysque mon cueur est tien.
		tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] += 1;
		indiv[k].fitness=fitness(k);
		if (maxfitness<indiv[k].fitness) maxfitness=indiv[k].fitness;
		}
	}

	  //replace by crossover

for(truc=0;truc<PCROSS;truc++)
     {
      i=(POPULATION-1)*randfloat0();
      j=(POPULATION-1)*randfloat0();
      if(j!=i)
        {
		if(indiv[i].fitness /(.01+tabsharing[(int)(indiv[i].xc/SHARINGRADIUS)][(int)(indiv[i].yc/SHARINGRADIUS)])
			<indiv[j].fitness /(.01+tabsharing[(int)(indiv[j].xc/SHARINGRADIUS)][(int)(indiv[j].yc/SHARINGRADIUS)]))
		{k=i;l=j;} else {k=j;l=i;}
			//we call k the loser and l the winner
		tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] -= 1;
			i=l; j=(POPULATION-1)*randfloat0();
//the winner l becomes the parent i, the other parent j is drawn randomly
			a=0;b=0;c=0;
         while ((c<MINRADIUS)||(c>MAXRADIUS)||(a<c)||(a>width-c)||(b<c)||(b>height-c))
		 {lambda=.5+randfloat1();                 // choose lambda between -0.5 and +1.5
		  a=(lambda*indiv[i].xc+(1-lambda)*indiv[j].xc);
          b=(lambda*indiv[i].yc+(1-lambda)*indiv[j].yc);
		  c=(lambda*indiv[i].r +(1-lambda)*indiv[j].r );
		 }
		indiv[k].xc=a;
		indiv[k].yc=b;
		indiv[k].r=c;
        indiv[k].fitness=fitness(k);
		tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS] +=1;
		if (maxfitness<indiv[k].fitness) maxfitness=indiv[k].fitness;
     	}
     }

if (generation==10000*(int)(generation/10000)) printf("\ngeneration no %d\tmaxfitness=%f",generation, maxfitness);
  }

/* now we have a population at the end of generations, what shall we do with it? */

for (i=1;i<width-1;i++) for (j=1;j<height-1;j++) ioutput[i][j]=iinput[i][j];

//exhibit the best solutions at end of algorithm

//printf ("\n evolution terminated, this is the best solutions I found:");
averagefitness=0;
for (k=0;k<POPULATION;k++) averagefitness += (indiv[k].fitness)*(indiv[k].fitness);
  
averagefitness = (averagefitness/POPULATION);
	
coefficient=.95; coeff=1.; //mistake was here in version 23
counter=0;
while ((counter<numberofcircles)&&(coeff>0.01))
	{
	 coeff *= coefficient;
  
  // MODIFICATION
  
   //  printf("counter = %d  coeff= %f\n",counter, coeff);
     
     
     
	 for(k=0;k<POPULATION;k++)
	 if((indiv[k].fitness>coeff*maxfitness)
		&&(indiv[k].retained==0)
		&&(tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS]>0)
	   )
	   {indiv[k].retained=1;
		counter++;
		tabsharing[indiv[k].xc/SHARINGRADIUS][indiv[k].yc/SHARINGRADIUS]=0;
	   }
	}
for(k=0;k<POPULATION;k++) if (indiv[k].retained==1)
{
  printf("\nindiv no. %d\tfitness=%f\txc=%d\tyc=%d\tr=%d",k,indiv[k].fitness,indiv[k].xc,indiv[k].yc,indiv[k].r);
  printf("%d, ",indiv[k].r);
  radius[k] = indiv[k].r;
  stat[radius[k]] ++;
for(i=1;i<width-1;i++) for(j=1;j<height-1;j++)
  {a=indiv[k].xc;b=indiv[k].yc;c=indiv[k].r;d=(a-i)*(a-i)+(b-j)*(b-j)-c*c;
   if(d-c<2) if(d+c>2)ioutput[i][j]=255;
   //if(d*d<4)ioutput[i][j]=255;
  }
}
fp = fopen("statistic.txt","wb");
for(k=0;k<POPULATION;k++) if (indiv[k].retained==1)fprintf(fp, "%d, ",radius[k]);
fclose(fp);


//MODIFICATION

//for(k=MINRADIUS;k<=MAXRADIUS;k++) printf ("\n radius = %d # of circles %d", k, stat[k]);

printf("\nTerminated: %d evaluations\n",nb_evaluations);
}
/*******************************************/
void showcircles()
{
}
/*******************************************/
float fitness(int i)
{

float sum, squaresum, standdev;
    float inside=0., squareinside=0., outside=0;;
int x1,y1,x2,y2,x3,y3;
sum=0.; squaresum=0.; standdev=0.; c=0.;
nb_evaluations++;                                   //just a counter
x=indiv[i].xc; y=indiv[i].yc; r=indiv[i].r;
if ((x<=r+5)||(x>=width-r-6)||(y<=r+4)||(y>=height-r-6)) c=0.;  // centre too close to the edges
else
  {
	c=2*SAMPLES*sqrt(r);
    
	for(j=0;j<SAMPLES;j++)
	{
// -1- bright edge
      theta=(int) (627.*randfloat0());
      x1=(int)((float)x+r*COS[theta]);
	  y1=(int)((float)y+r*SIN[theta]); //this is a random point on the circle
	  c += edge * (iinput[x1][y1]); // high grey level on edge contributes positively to fitness
      
// -2- high border contrast

      c += border * igrad[x1][y1]; // bonus if high contrast on edge of circle
    
// -3- contrast between inside and edge
      a=randfloat0(); /*a=sqrt(1-a*a);*/ a=0.66*r*a;   //a is between 0 and r
	  x1=(int)(x+a*COS[theta]);y1=(int)(y+a*SIN[theta]); //inside  disc
	  x2=(int)(x+(r)*COS[theta]);y2=(int)(y+(r)*SIN[theta]);  //on edge of disc along same radial
        inside = iinput[x1][y1];  //inside
        outside = iinput[x2][y2];  //on edge (not outside)
      c += contrast * (outside - inside);  // happy if outside brighter than inside

// -4- homogeneity between two points inside the disc, at different angles and -un-equal distance from the centre
      theta=(int) (627.*randfloat0());   //a new theta
      a *= .8 * randfloat0() * r ;
	  x3=(int)(x+a*COS[theta]);y3=(int)(y+a*SIN[theta]);  //second point inside
      c -= homogeneity * (iinput[x1][y1]-iinput[x3][y3])*(iinput[x1][y1]-iinput[x3][y3]);
           // terribly angry if inside points have very different brightness

// -5- (bonus if inside is dark) penalty if inside is bright
      c -= darkness*iinput[x1][y1];
      
// -6- bonus if diameter is close to 5 or close to 10
     //c += prefdiameter / (1+(r-5)*(r-5));
      c += prefdiameter / (1+(r-100)*(r-100));
    }
  }
// return c*c; //in order to reduce influence of sharing
return c;

}

/*******************************************/
float randfloat0(void)
{
  float alea, numerator, denominator;

  numerator = rand() & RAND_MAX;
  denominator = RAND_MAX;

  alea = numerator/denominator;
  return alea;  /* equiprobable over [0,1] */
}
/**************************************************************** RANDFLOAT1 */
float randfloat1(void)
{
  float alea, numerator, denominator;

  numerator = rand() & RAND_MAX;
  denominator = RAND_MAX;

  alea = 2. *(numerator/denominator) - 1.;
  alea = alea * (.2 + .8 * alea * alea ) ;
  return alea;  /* Gaussian-like over [-1,1], peak on 0 */
}

/*******************************************/

int main()
{
readinput();
initialise();
process();
writeoutput();
}
