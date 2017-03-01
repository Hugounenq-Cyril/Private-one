#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 140

void print(float a[N])
{
    printf("(");
    int i;
    for(i=0;i<N-1;i++)
    {
        printf("%f,",a[i]);
    }
    printf("%f)\n",a[N-1]);
}

float abs(float t){

	if(t>=0) return t;
	else return -t;

}
int countnegs(float a[])
{
    int negs = 0;
    int i;
    for(i=0;i<N;i++)
    {
        if(a[i]<0) negs++;
    }
    return negs;
}


void convex(float t, float u[], float res[])
{
	//f(t)=(1-|t|)u+tv
	int i;
	for(i=0;i<N;i++)
	{
		res[i]=(1-abs(t))*u[i]+t*1;
	}

}

float t0rec(float t0,float t1, float u[])
{
    float aux[N];
    float param = (t0+t1)/2; 
    
    convex(param,u,aux);
    
    printf("%d\n",countnegs(aux));
    if (countnegs(aux)==0){ return t0rec(t0,param,u);}

    if(countnegs(aux)==1){return param;}
    
    else { return t0rec(param,t1,u);}

    

return 55;

}

main() {
    
    float ar[N];    
    int i;
    float S=0;
    int j;
    
for(j = 0;j<10;j++)
{

    for(i=0;i<N;i++)
    {
        ar[i]=rand()%1000-500;
        S+=ar[i];
    }
    for(i=0;i<N;i++)
    {
        ar[i]-=S/N;
    }
    
    
    
    
    
    float Norm = 0;
    for(i=0;i<N;i++)
    {
        Norm+=ar[i]*ar[i];
    }
    Norm=sqrt(Norm);
    
    for(i=0;i<N;i++)
    {
        ar[i]/=Norm;
        ar[i]*=sqrt(N);
    }
    
    //print(ar);
    S=0;
    Norm=0;
    for(i=0;i<N;i++)
    {
        S+=ar[i];
        Norm+=ar[i]*ar[i];
    }
    printf("S=%f N=%f\n",S,Norm);

    float f[N];
    
    float beginning = 0.6488;
    float end = 0.65;
    float largo = end-beginning;
    float iterations = 10;
    float step = largo/iterations;
    float t;

    for(i=0;i<=iterations;i++)
	{
	    t=beginning+i*step;
	    convex(t,ar,f);
	    //print(f);
	    //printf("t = %f, negatives = %d\n",t,countnegs(f));	
    }

    printf("param = %f\n",t0rec(0.5,1,ar));
}

}