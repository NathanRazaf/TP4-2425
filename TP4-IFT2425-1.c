//------------------------------------------------------
// module  : Tp4-IFT2425-1.c
// author  : 
// date    : 
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//  

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <new>
#include <unistd.h>

/************************************************************************/
/* WINDOWS						          	*/
/************************************************************************/
#include <X11/Xutil.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;

//------------------------------------------------
// DEFINITIONS -----------------------------------                       
//------------------------------------------------
#define CARRE(X) ((X)*(X))

#define OUTPUT_FILE "Tp4-Img-I.pgm"
#define VIEW_PGM    "xv" 
#define DEBUG 0

//-Cst-Modele
#define X_1 0.0
#define Y_1 1.0
#define X_2 -1.0/sqrt(2.0)
#define Y_2 -1.0/2.0
#define X_3 +1.0/2*sqrt(2.0)
#define Y_3 -1.0/2.0
#define C 0.25 
#define R 0.1
#define D 0.3

#define X_1_INI 0.2            
#define X_2_INI 0.0       
#define X_3_INI -1.6                  
#define X_4_INI 0.0 
 
//-Cst-Runge-Kutta
#define H            0.0001  
#define T_0          0.0                 
#define T_F         30.0 
#define NB_INTERV (T_F-T_0)/H
   
 //-Cst-Image                             
#define WIDTH  512
#define HEIGHT 512                  
#define MAX_X  4.0                
#define MAX_Y  4.0  
#define EVOL_GRAPH 3000
              
#define WHITE     255
#define GREYWHITE 230
#define GREY      200
#define GREYDARK  120
#define BLACK       0   

//------------------------------------------------
// GLOBAL CST ------------------------------------                       
//------------------------------------------------
float Xmin=0.0;
float Xmax=0.0;
float Ymin=0.0;
float Ymax=0.0; 

float xx_1=((WIDTH/MAX_X)*X_1)+(WIDTH/2);
float yy_1=(-(HEIGHT/MAX_Y)*Y_1)+(HEIGHT/2);
float xx_2=((WIDTH/MAX_X)*X_2)+(WIDTH/2);
float yy_2=(-(HEIGHT/MAX_Y)*Y_2)+(HEIGHT/2);
float xx_3=((WIDTH/MAX_X)*X_3)+(WIDTH/2);
float yy_3=(-(HEIGHT/MAX_Y)*Y_3)+(HEIGHT/2);

/************************************************************************/
/* OPEN_DISPLAY()							*/
/************************************************************************/
int open_display()
{
  if ((display=XOpenDisplay(NULL))==NULL)
   { printf("Connection impossible\n");
     return(-1); }

  else
   { screen_num=DefaultScreen(display);
     visual=DefaultVisual(display,screen_num);
     depth=DefaultDepth(display,screen_num);
     root=RootWindow(display,screen_num);
     return 0; }
}

/************************************************************************/
/* FABRIQUE_WINDOW()							*/
/* Cette fonction crée une fenetre X et l'affiche à l'écran.	        */
/************************************************************************/
Window fabrique_window(char *nom_fen,int x,int y,int width,int height,int zoom)
{
  Window                 win;
  XSizeHints      size_hints;
  XWMHints          wm_hints;
  XClassHint     class_hints;
  XTextProperty  windowName, iconName;

  char *name=nom_fen;

  if(zoom<0) { width/=-zoom; height/=-zoom; }
  if(zoom>0) { width*=zoom;  height*=zoom;  }

  win=XCreateSimpleWindow(display,root,x,y,width,height,1,0,255);

  size_hints.flags=PPosition|PSize|PMinSize;
  size_hints.min_width=width;
  size_hints.min_height=height;

  XStringListToTextProperty(&name,1,&windowName);
  XStringListToTextProperty(&name,1,&iconName);
  wm_hints.initial_state=NormalState;
  wm_hints.input=True;
  wm_hints.flags=StateHint|InputHint;
  class_hints.res_name=nom_fen;
  class_hints.res_class=nom_fen;

  XSetWMProperties(display,win,&windowName,&iconName,
                   NULL,0,&size_hints,&wm_hints,&class_hints);

  gc=XCreateGC(display,win,0,NULL);

  XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask| 
               ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask| 
               StructureNotifyMask);

  XMapWindow(display,win);
  return(win);
}

/****************************************************************************/
/* CREE_XIMAGE()							    */
/* Crée une XImage à partir d'un tableau de float                          */
/* L'image peut subir un zoom.						    */
/****************************************************************************/
XImage* cree_Ximage(float** mat,int z,int length,int width)
{
  int lgth,wdth,lig,col,zoom_col,zoom_lig;
  float somme;
  unsigned char	 pix;
  unsigned char* dat;
  XImage* imageX;

  /*Zoom positiv*/
  /*------------*/
  if (z>0)
  {
   lgth=length*z;
   wdth=width*z;

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
   { 
    pix=(unsigned char)mat[lig/z][col/z];
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
      { 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pix; 
       }
    }
  } /*--------------------------------------------------------*/

  /*Zoom negatifv*/
  /*------------*/
  else
  {
   z=-z;
   lgth=(length/z);
   wdth=(width/z);

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
   {  
    somme=0.0;
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
     somme+=mat[lig+zoom_lig][col+zoom_col];
           
     somme/=(z*z);    
     dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)somme; 
   }
  } /*--------------------------------------------------------*/

  imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
  return (imageX);
}

//------------------------------------------------
// FUNCTIONS -------------------------------------                       
//------------------------------------------------
//-------------------------//
//-- Matrice de Double ----//
//-------------------------//
//---------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* dmatrix_allocate_1d(int hsize)
 {
  float* matrix;
  matrix=new float[hsize]; return matrix; }

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** dmatrix_allocate_2d(int vsize,int hsize)
 {
  float** matrix;
  float *imptr;

  matrix=new float*[vsize];
  imptr=new float[(hsize)*(vsize)];
  for(int i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_dmatrix_1d(float* pmat)
{ delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_dmatrix_2d(float** pmat)
{ delete[] (pmat[0]);
  delete[] pmat;}

//----------------------------------------------------------
// SaveImagePgm                       
//----------------------------------------------------------                
void SaveImagePgm(char* name,float** mat,int lgth,int wdth)
{
 int i,j;
 char buff[300];
 FILE* fic;

  //--extension--
  strcpy(buff,name);

  //--ouverture fichier--
  fic=fopen(buff,"wb");
    if (fic==NULL) 
        { printf("Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format pgm\n",buff);

  //--sauvegarde de l'entete--
  fprintf(fic,"P5");
  fprintf(fic,"\n# IMG Module");
  fprintf(fic,"\n%d %d",wdth,lgth);
  fprintf(fic,"\n255\n");

  //--enregistrement--
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
	fprintf(fic,"%c",(char)mat[i][j]);
   
  //--fermeture fichier--
  fclose(fic); 
}

//------------------------------------------------------------------------
// plot_point
//
// Affiche entre x dans [-MAX_X/2  MAX_X/2]
//               y dans [-MAX_Y/2  MAX_Y/2]                
//------------------------------------------------------------------------
void plot_point(float** MatPts,float** MatPict,int NbPts)
{
 int x_co,y_co;
 int i,j,k;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)  MatPict[i][j]=GREYWHITE;

 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 for(k=0;k<NbPts;k++)
    { x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if (DEBUG) printf("[%d::%d]",x_co,y_co); 
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
	 MatPict[y_co][x_co]=BLACK; 
    }
}

//------------------------------------------------------------------------
// Fill_Pict
//------------------------------------------------------------------------
void Fill_Pict(float** MatPts,float** MatPict,int PtsNumber,int NbPts)
{
 int i,j;
 int x_co,y_co;
 int k,k_Init,k_End;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if (MatPict[i][j]!=GREYWHITE) MatPict[i][j]=GREY;
     if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 k_Init=PtsNumber;
 k_End=(k_Init+EVOL_GRAPH)%NbPts;
 for(k=k_Init;k<k_End;k++)
    { k=(k%NbPts);
      x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
         MatPict[y_co][x_co]=BLACK; }
}


//------------------------------------------------
// FONCTIONS TPs----------------------------------                      
//------------------------------------------------
      
// Fonction pour calculer les dérivées du système d'équations différentielles
void compute_derivatives(float t, float x, float y, float vx, float vy, float* dvx, float* dvy)
{
    float sum_x = 0;
    float sum_y = 0;
    float dist1, dist2, dist3;

    // Calcul des distances au cube
    dist1 = pow(sqrt(CARRE(X_1 - x) + CARRE(Y_1 - y) + CARRE(D)), 3);
    dist2 = pow(sqrt(CARRE(X_2 - x) + CARRE(Y_2 - y) + CARRE(D)), 3);
    dist3 = pow(sqrt(CARRE(X_3 - x) + CARRE(Y_3 - y) + CARRE(D)), 3);

    // Somme des forces d'attraction des aimants
    sum_x = (X_1 - x) / dist1 + (X_2 - x) / dist2 + (X_3 - x) / dist3;
    sum_y = (Y_1 - y) / dist1 + (Y_2 - y) / dist2 + (Y_3 - y) / dist3;

    // Calcul des dérivées secondes (accélérations)
    *dvx = -R * vx + sum_x - C * x;
    *dvy = -R * vy + sum_y - C * y;
}

// Implémentation de la méthode de Runge-Kutta Fehlberg
void runge_kutta_fehlberg(float** MatPts, int nb_points)
{
    float t = T_0;
    float x = X_1_INI;    // Position initiale en x
    float y = X_3_INI;    // Position initiale en y
    float vx = X_2_INI;   // Vitesse initiale en x
    float vy = X_4_INI;   // Vitesse initiale en y

    float k1_x, k1_y, k1_vx, k1_vy;
    float k2_x, k2_y, k2_vx, k2_vy;
    float k3_x, k3_y, k3_vx, k3_vy;
    float k4_x, k4_y, k4_vx, k4_vy;
    float k5_x, k5_y, k5_vx, k5_vy;
    float k6_x, k6_y, k6_vx, k6_vy;

    float dvx, dvy;

    // Initialisation des statistiques
    Xmin = Xmax = x;
    Ymin = Ymax = y;

    // Stocker la position initiale
    MatPts[0][0] = x;
    MatPts[0][1] = y;

    // Boucle principale de la méthode RKF
    for (int i = 1; i < nb_points; i++)
    {
        // Calcul de k1
        k1_x = H * vx;
        k1_y = H * vy;
        compute_derivatives(t, x, y, vx, vy, &dvx, &dvy);
        k1_vx = H * dvx;
        k1_vy = H * dvy;

        // Calcul de k2
        k2_x = H * (vx + k1_vx/4);
        k2_y = H * (vy + k1_vy/4);
        compute_derivatives(t + H/4, x + k1_x/4, y + k1_y/4, vx + k1_vx/4, vy + k1_vy/4, &dvx, &dvy);
        k2_vx = H * dvx;
        k2_vy = H * dvy;

        // Calcul de k3
        k3_x = H * (vx + 3*k1_vx/32 + 9*k2_vx/32);
        k3_y = H * (vy + 3*k1_vy/32 + 9*k2_vy/32);
        compute_derivatives(t + 3*H/8, x + 3*k1_x/32 + 9*k2_x/32, y + 3*k1_y/32 + 9*k2_y/32,
                         vx + 3*k1_vx/32 + 9*k2_vx/32, vy + 3*k1_vy/32 + 9*k2_vy/32, &dvx, &dvy);
        k3_vx = H * dvx;
        k3_vy = H * dvy;

        // Calcul de k4
        k4_x = H * (vx + 1932*k1_vx/2197 - 7200*k2_vx/2197 + 7296*k3_vx/2197);
        k4_y = H * (vy + 1932*k1_vy/2197 - 7200*k2_vy/2197 + 7296*k3_vy/2197);
        compute_derivatives(t + 12*H/13,
                         x + 1932*k1_x/2197 - 7200*k2_x/2197 + 7296*k3_x/2197,
                         y + 1932*k1_y/2197 - 7200*k2_y/2197 + 7296*k3_y/2197,
                         vx + 1932*k1_vx/2197 - 7200*k2_vx/2197 + 7296*k3_vx/2197,
                         vy + 1932*k1_vy/2197 - 7200*k2_vy/2197 + 7296*k3_vy/2197, &dvx, &dvy);
        k4_vx = H * dvx;
        k4_vy = H * dvy;

        // Calcul de k5
        k5_x = H * (vx + 439*k1_vx/216 - 8*k2_vx + 3680*k3_vx/513 - 845*k4_vx/4104);
        k5_y = H * (vy + 439*k1_vy/216 - 8*k2_vy + 3680*k3_vy/513 - 845*k4_vy/4104);
        compute_derivatives(t + H,
                         x + 439*k1_x/216 - 8*k2_x + 3680*k3_x/513 - 845*k4_x/4104,
                         y + 439*k1_y/216 - 8*k2_y + 3680*k3_y/513 - 845*k4_y/4104,
                         vx + 439*k1_vx/216 - 8*k2_vx + 3680*k3_vx/513 - 845*k4_vx/4104,
                         vy + 439*k1_vy/216 - 8*k2_vy + 3680*k3_vy/513 - 845*k4_vy/4104, &dvx, &dvy);
        k5_vx = H * dvx;
        k5_vy = H * dvy;

        // Calcul de k6
        k6_x = H * (vx - 8*k1_vx/27 + 2*k2_vx - 3544*k3_vx/2565 + 1859*k4_vx/4104 - 11*k5_vx/40);
        k6_y = H * (vy - 8*k1_vy/27 + 2*k2_vy - 3544*k3_vy/2565 + 1859*k4_vy/4104 - 11*k5_vy/40);
        compute_derivatives(t + H/2,
                         x - 8*k1_x/27 + 2*k2_x - 3544*k3_x/2565 + 1859*k4_x/4104 - 11*k5_x/40,
                         y - 8*k1_y/27 + 2*k2_y - 3544*k3_y/2565 + 1859*k4_y/4104 - 11*k5_y/40,
                         vx - 8*k1_vx/27 + 2*k2_vx - 3544*k3_vx/2565 + 1859*k4_vx/4104 - 11*k5_vx/40,
                         vy - 8*k1_vy/27 + 2*k2_vy - 3544*k3_vy/2565 + 1859*k4_vy/4104 - 11*k5_vy/40, &dvx, &dvy);
        k6_vx = H * dvx;
        k6_vy = H * dvy;

        // Mise à jour des positions et vitesses
        x = x + (16*k1_x/135 + 6656*k3_x/12825 + 28561*k4_x/56430 - 9*k5_x/50 + 2*k6_x/55);
        y = y + (16*k1_y/135 + 6656*k3_y/12825 + 28561*k4_y/56430 - 9*k5_y/50 + 2*k6_y/55);
        vx = vx + (16*k1_vx/135 + 6656*k3_vx/12825 + 28561*k4_vx/56430 - 9*k5_vx/50 + 2*k6_vx/55);
        vy = vy + (16*k1_vy/135 + 6656*k3_vy/12825 + 28561*k4_vy/56430 - 9*k5_vy/50 + 2*k6_vy/55);

        // Avancer le temps
        t += H;

        // Stocker les positions
        MatPts[i][0] = x;
        MatPts[i][1] = y;

        // Mettre à jour les statistiques
        if (x < Xmin) Xmin = x;
        if (x > Xmax) Xmax = x;
        if (y < Ymin) Ymin = y;
        if (y > Ymax) Ymax = y;
    }
}

//----------------------------------------------------------
//----------------------------------------------------------
// MAIN  
//----------------------------------------------------------
//----------------------------------------------------------
int main (int argc, char **argv)
{
  int i,j,k;
  int flag_graph;
  int zoom;

  XEvent ev;
  Window win_ppicture;
  XImage *x_ppicture;
  char   nomfen_ppicture[100]; 
  char BufSystVisu[100];

  //>AllocMemory
  float** MatPict=dmatrix_allocate_2d(HEIGHT,WIDTH);
  float** MatPts=dmatrix_allocate_2d((int)(NB_INTERV),2);
  
  //>Init
  for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
  for(i=0;i<2;i++) for(j=0;j<(int)(NB_INTERV);j++) MatPts[i][j]=0.0;
  flag_graph=1;
  zoom=1;


  //---------------------------------------------------------------------
  //>Question 1 
  //---------------------------------------------------------------------  

  //Il faut travailler ici ...et dans > // FONCTIONS TPs

  //Un exemple ou la matrice de points est remplie
  //par une courbe donné par l'équation d'en bas... et non pas par 
  //la solution de l'équation différentielle
 
  // Appel de la fonction Runge-Kutta Fehlberg pour résoudre l'équation différentielle
  runge_kutta_fehlberg(MatPts, (int)(NB_INTERV));

  //--Fin Question 1-----------------------------------------------------


  //>Affichage des Points dans MatPict
  plot_point(MatPts,MatPict,(int)(NB_INTERV));

  //>Save&Visu de MatPict
  SaveImagePgm((char*)OUTPUT_FILE,MatPict,HEIGHT,WIDTH);
  strcpy(BufSystVisu,VIEW_PGM);
  strcat(BufSystVisu," "); 	
  strcat(BufSystVisu,OUTPUT_FILE);
  strcat(BufSystVisu," &"); 
  system(BufSystVisu);

  //>Affiche Statistique
  printf("\n\n Stat:  Xmin=[%.2f] Xmax=[%.2f] Ymin=[%.2f] Ymax=[%.2f]\n",Xmin,Xmax,Ymin,Ymax);


 //--------------------------------------------------------------------------------
 //-------------- visu sous XWINDOW de l'évolution de MatPts ----------------------
 //--------------------------------------------------------------------------------
 if (flag_graph)
 {
 //>Uuverture Session Graphique
 if (open_display()<0) printf(" Impossible d'ouvrir une session graphique");
 sprintf(nomfen_ppicture,"Évolution du Graphe");
 win_ppicture=fabrique_window(nomfen_ppicture,10,10,HEIGHT,WIDTH,zoom);
 x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);

 printf("\n\n Pour quitter,appuyer sur la barre d'espace");
 fflush(stdout);

 //>Boucle Evolution
  for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
  for(k=0;;)
     {   
       k=((k+EVOL_GRAPH)%(int)(NB_INTERV));
       Fill_Pict(MatPts,MatPict,k,(int)(NB_INTERV));
       XDestroyImage(x_ppicture);
       x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);
       XPutImage(display,win_ppicture,gc,x_ppicture,0,0,0,0,x_ppicture->width,x_ppicture->height); 
       usleep(10000);  //si votre machine est lente mettre un nombre plus petit
     }
 } 
       
 //>Retour  
 printf("\n Fini... \n\n\n");
 return 0;
}




