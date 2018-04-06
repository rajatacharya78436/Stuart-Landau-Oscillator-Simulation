/* Stuart-landau oscillator ***********/

#include<stdio.h>

#include<math.h>

#include<stdlib.h>

#include<limits.h>

#include<time.h>

#define pi 4.0*atan(1.0)

#define RANDMAX 2147483647

void rk41(float t,float x[],float h);

float F1(float t,float x[],int i);

int NSYS=3;

int n=3;

float dt = 0.01;

float e,g,et,om,a,sig,rm;

int j,ip1,ip2;


void main()

{ 

    
FILE *a_out0,*a_out1;

     a_out0=fopen("sl_t","w");
 
    a_out1=fopen("sl_s","w");

   
  int i,k,na,l,ip,is,iq,ic,NP,im,ij,l1,j1,k1,k2,k3;

     na=150000;

   
  float t = 0.0,x[n];

    
 float max=1,min=-1; 
 
  for(i=1;i<=NSYS;i++)

   {
     
   
    x[i]=min+(max-min)*(double)(rand())/(double)(RANDMAX);
 
   }
 

   /***************************/


    sig=4;
 rm=1.35;

    e=((sig*(rm*rm))/2);

    g=(-sig/4);

    et=0.1;
 
    om=3;
 
    a=0.8;




    for(i=1;i<=na;i++)
  
    {

                       
     rk41(t,x,dt);

 
            if(i>=130000)
 
        
     {
         
         
		fprintf(a_out0,"%f\t%f\t%f\t%f\n",t,x[1],x[2],x[3]);

 
        
     }   

      
 
    t=t+dt;
   
 } 
        
   
     
 

  
   


}
          
         
 


 
/*********************-Equation------------------------*/

float F1(float t,float x[],int i)

  {
         
        
       
      if(i==1)
          return(x[3]*x[1]-om*x[2]+2*pow(x[1],3)+2*x[1]*pow(x[2],2)-(e*x[2]*pow(x[1],2))-(e*pow(x[2],3))-pow(x[1],5)-(x[1]*pow(x[2],4))-(2*pow(x[1],3))*pow(x[2],2)-(g*x[2]*pow(x[1],4))-(g*pow(x[2],5))-(2*g*pow(x[1],2)*pow(x[2],3)));
 

      if(i==2) 
          return(x[3]*x[2]+om*x[1]+(2*x[2]*pow(x[1],2))+2*pow(x[2],3)+(e*pow(x[1],3))+(e*x[1]*pow(x[2],2))-(x[2]*pow(x[1],4))-pow(x[2],5)-(2*pow(x[1],2)*pow(x[2],3))+g*pow(x[1],5)+g*x[1]*pow(x[2],4)+2*g*pow(x[1],3)*pow(x[2],2));   
   

      if(i==3)
          return(et*(a-(pow(x[1],2)+pow(x[2],2))));

  } 
  
 

   


/******************** --RK4 FUNCTION -------------------------------------------*/



void rk41(float t,float x[],float step)


 {
  
    float h = step/2.0;

    float t1[n],t2[n],t3[n],k1[n],k2[n],k3[n],k4[n];
    
    int i;
 
      for(i=1;i<=n;i++) 
            t1[i] = x[i] + 0.5*(k1[i] = step*F1(t,x,i));
 
      for(i=1;i<=n;i++)
            t2[i] = x[i] + 0.5*(k2[i] = step*F1(t+h,t1,i)); 

      for(i=1;i<=n;i++) 
            t3[i] = x[i] + (k3[i] = step*F1(t+h,t2,i));
 
      for(i=1;i<=n;i++) 
            k4[i] = dt*F1(t+step,t3,i);
  

      for(i=1;i<=n;i++)
            x[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
 } 
















