#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {

  int g,h,i,j,K,M,poss,maxIndex,counter,counterPointer,first,last,dummyInt; bool val;
  double* T, *meanDiff, *vec,*vec2, *dummyPointer, *dummyArray,*dummyArray2, *dummyArrayPattern, *vecPattern;
  double threshold, scale, shift, maxDist=-1,dist,adummy,bdummy,cdummy;

  bool better=1;

  // Test number of parameters.
  if (nrhs != 3 || nlhs != 3) {
    mexWarnMsgTxt("Usage: [Cand,Poss, dist] = getVertices(T',meanDiff, V)\n");
    return;
  }

  // Parse parameters
  K   = (int)mxGetM(prhs[0]); // number of rows of T
  M   = (int)mxGetN(prhs[0]); // number of cols of T
  T = mxGetPr(prhs[0]);

  if( ((int)mxGetM(prhs[1])) != M)
   mexErrMsgTxt("Number of columns of meanDiff and columns of TT are not the same");

  meanDiff = mxGetPr(prhs[1]);
  
  double* V = mxGetPr(prhs[2]);
  threshold = 0.5 * (V[0] + V[1]);
  scale = V[1] - V[0];
  shift = V[0];
  
  //mexPrintf("threhold: %f\n",threshold);
  //mexPrintf("scale: %f\n",scale);
   //mexPrintf("shift: %f\n",shift);
  
   // Allocate memory for output
  plhs[0] = mxCreateDoubleMatrix(M, K+1, mxREAL);
  double* Vertices = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(K, K+1, mxREAL);
  double* VerticesPattern = mxGetPr(plhs[1]);
  
  plhs[2] = mxCreateDoubleMatrix(1, K+1, mxREAL);
  double* distances_out = mxGetPr(plhs[2]);

  if(Vertices==NULL || VerticesPattern==NULL) mexErrMsgTxt("Pointer of Vertices failed");

  double* ArrayPatternGlobal = new double[(K+2)*K];
  double** ArrayPattern = new double*[K+2];
  for(i=0;i<K+2;i++)
   { ArrayPattern[i]=ArrayPatternGlobal+i*K;  }
  dummyArrayPattern=ArrayPattern[K+1];

  double* ArrayGlobal=new double[(K+3)*M];
  if(ArrayGlobal==NULL) mexErrMsgTxt("Pointer of ArrayGlobal failed");
  double** Array =new double*[K+3];
  for(i=0;i<K+3;i++)
   { Array[i]=ArrayGlobal+i*M; if(Array[i]==NULL) mexErrMsgTxt("Pointer of array failed");}

  dummyArray =Array[K+1];
  dummyArray2=Array[K+2];
  if(dummyArray==NULL) mexErrMsgTxt("Pointer of dummyArray failed");

  double* distances=new double[K+1];//(double*)mxMalloc((K+1)*sizeof(double)); 
  //
  if(distances==NULL) mexErrMsgTxt("Pointer of distances failed");

  bool bit;
  int* indices=new int[K+1];
  int* diff=new int[K+1];

  int base;
  if(K<=10)
   base=5;
  else
   base=10;
  if(M>5000)
   base=8;
  if(M>100000)
   base=6;
  if(base>K)
   base=K-1;
  int powerbase=1;
  for(i=0;i<base;i++)
   powerbase*=2;
  //mexPrintf("Powerbase: %i\n",powerbase);


  double* PreCompVals=new double[powerbase*M];
  if(PreCompVals==NULL) mexErrMsgTxt("Pointer of precomputed values failed");

  // number of possibilities
  poss=1;
  for(i=0;i<K;i++)
   poss*=2;

  
  /*for(i=0; i<K+1; i++)
  {
	dummyArray=Array[i];
    for(j=0; j<M; j++)
	{
	  dummyArray[j]=meanDiff[j];
	}
  }*/

  clock_t firstpart, secondpart,total_start;
  total_start = clock();
  
  for(i=0; i<(K+1); i++) 
  //for(i=poss-1; i>poss-1-(K+1); i--)
  {
	vec=Array[i]; dist=0; counter=0; vecPattern=ArrayPattern[i];
	//vec=Array[i-(poss-1-(K+1))-1]; dist=0;
	for(h=0; h<K; h++)
    {  
	  bit=(i >> h) & 1; // get the value of the h-bit of i
	  //if(bits[h] & minBit==-1)
	  // minBit=h; 
	  if(bit)
	   { indices[counter]=h; counter++; vecPattern[h]=1; }
	  else
	   { vecPattern[h]=0; }
	}
	for(j=0; j<M; j++)
	{
	  dummyPointer=T+j*K;  // go to the j-th column of TT
	  adummy=meanDiff[j];
	  for(h=0;h<counter;h++)
	   adummy+=dummyPointer[indices[h]];

	  /*for(h=minBit; h<maxBit; h++)
      {  
	    //val=(i >> h) & 1; // get the value of the h-bit of 
		if(j==0){ mexPrintf("Bit: %i, Value: %i\n",h,val); }
	    if(bits[h])           // if yes add the vector     
	     adummy+=dummyPointer[h];
	  }*/
      // threshold and compute distance
	  vec[j]=adummy;
	  if(adummy>threshold)
	   { bdummy=scale + shift - adummy; dist+=bdummy*bdummy; }
	  else
	   { bdummy = shift - adummy; dist+=bdummy*bdummy; }
	}
	distances[i]=dist;
	//distances[i-(poss-1-(K+1))-1]=dist;
	//mexPrintf("Iteration %i: dist: %1.6f\n",i,dist);
	if(dist>maxDist)
     maxDist=dist;
  }
  //mexPrintf("Iteration %i: Initial maxdist: %1.6f, - poss: %i\n",i,maxDist,poss
  

  for(i=(K+1); i<powerbase; i++)
  {
    vec=dummyArray;  dist=0; better=true; counter=0;//minBit=0; maxBit=0;
	for(h=0; h<K; h++)
    {  
	  bit=(i >> h) & 1; // get the value of the h-bit of i
	  if(bit)
	   { if(counter>0){ diff[counter-1]=h-last; last=h;}else{first=h; last=h;} counter++;}
	}
    diff[counter-1]=K+first-last; //+indices[0]-indices[counter-1];
	//indices[counter]=K+indices[0];
	//for(h=0;h<counter;h++)
    // diff[h]=indices[h+1]-indices[h];

	dummyPointer=T+first;
	for(j=0; j<M; j++)
	{
      adummy=meanDiff[j];
	  //counterPointer;
	  //dummyPointer=T+j*K;
	  for(h=0;h<counter;h++)
	  { adummy+=(*dummyPointer); dummyPointer+=diff[h]; }//dummyPointer[indices[h]];
	  //counterPointer+=M;
	  /*for(h=minBit; h<maxBit; h++)
      {  
	    //val=(i >> h) & 1; // get the value of the h-bit of i
	    if(bits[h])           // if yes add the vector     
	     adummy+=dummyPointer[h];
	  }*/
      // threshold and compute distance
	  if(adummy>threshold)
	   { bdummy=scale + shift - adummy; dist+=bdummy*bdummy; }
	  else
	   { bdummy = shift - adummy; dist+=bdummy*bdummy; }
	  if(dist>maxDist )
	   {better=false; break; }
	  vec[j]=adummy;
	}
	if(better)
	{
      vecPattern=dummyArrayPattern;
      for(h=0; h<K; h++)
      {  
	    bit=(i >> h) & 1; // get the value of the h-bit of i
	    if(bit)
	     { vecPattern[h]=1; }
	    else
	     { vecPattern[h]=0; }
	  }
	  for(j=0; j<K+1; j++)
	  {
		if(distances[j]==maxDist)
		{ distances[j]=dist; dummyPointer=Array[j]; Array[j]=dummyArray; dummyArray=dummyPointer; 
		  dummyPointer=ArrayPattern[j]; ArrayPattern[j]=dummyArrayPattern; dummyArrayPattern=dummyPointer; break;}
	  }
	  maxDist=-1;
	  for(j=0; j<K+1; j++)
	  {
		if(distances[j]>maxDist)
		{ maxDist=distances[j]; }
	  }
	  //mexPrintf("Iteration: %i, New max dist: %1.6f\n",i,maxDist);
	}
  }
  firstpart=clock()-total_start;

  
  counter=0; counterPointer=0;
  for(i=0; i<powerbase; i++)
  {
	counter=0; dummyInt=i;
	for(h=0; h<base; h++)
    {  
	  //bit=(i >> h) & 1; // get the value of the h-bit of i
	  if(dummyInt & 1)
	   { if(counter>0){ diff[counter-1]=h-last; last=h;}else{first=h; last=h;} counter++;}
	  dummyInt=dummyInt >> 1;
	}
	//mexPrintf("counter: %i, %i\n",i,counter);
	if(counter>0)
	{
	  diff[counter-1]=K+first-last; 
	  dummyPointer=T+first;
	  counterPointer=i*M;
	  for(j=0; j<M; j++)
      {
        adummy=meanDiff[j];
	    for(h=0;h<counter;h++)
	      { adummy+=(*dummyPointer); dummyPointer+=diff[h]; }//dummyPointer[indices[h]];
	    PreCompVals[counterPointer]=adummy;
		//if(j==M-1) mexPrintf("precomp: %1.6f, %i\n",adummy,counterPointer);
	    counterPointer++;
	  }
	}
	else
	{
	  counterPointer=i*M;
	  for(j=0; j<M; j++)
	   { PreCompVals[counterPointer]=meanDiff[j]; counterPointer++;}
	}
  }

  for(i=powerbase; i<poss; i+=powerbase)
  {
    vec=dummyArray2; better=true; counter=0;//minBit=0; maxBit=0;
	/*dummyInt= (i>> base);
	for(h=base; h<K; h++)
    {  
	  //bit=(i >> h) & 1; // get the value of the h-bit of i
	  if(dummyInt & 1)
	   { if(counter>0){ diff[counter-1]=h-last; last=h;}else{first=h; last=h;} counter++;}
	  dummyInt=dummyInt >> 1;
	}
	//if(counter>0)
     diff[counter-1]=K+first-last; //+indices[0]-indices[counter-1];*/
	//else
	// first=0;
	//indices[counter]=K+indices[0];
	//for(h=0;h<counter;h++)
    // diff[h]=indices[h+1]-indices[h];
	//dummyPointer=T+first;
	dummyPointer=T+base; first= (i >> base); 
	for(j=0; j<M; j++)
	{
      adummy=0; //meanDiff[j];
	  //counterPointer;
	  //dummyPointer=T+j*K;
	  dummyInt= first; 
	  for(h=base; h<K; h++)
	  {
		if(dummyInt & 1) adummy+=(*dummyPointer);
		dummyPointer++; dummyInt=dummyInt >> 1;
	  }
	  dummyPointer+=base;
	  //for(h=0;h<counter;h++)
	  //  { adummy+=*dummyPointer; dummyPointer+=diff[h]; }//dummyPointer[indices[h]];
	  vec[j]=adummy;
	}
	counter=0; 
	for(g=0; g<powerbase; g++)
    {
	   dist=0; vec2=dummyArray; counter=g*M; better=true;
	   for(j=0; j<M; j++)
	   {
		 adummy=vec[j]+PreCompVals[counter]; 
	     if(adummy>threshold)
	      { bdummy=scale + shift -adummy; dist+=bdummy*bdummy; }
	     else
	      { bdummy = shift - adummy; dist+=bdummy*bdummy; }
	     if(dist>maxDist )
	      {better=false; break; }
		 vec2[j]=adummy;
		 counter++;
	   }
	   if(better)
	   {
         vecPattern=dummyArrayPattern;
         for(h=0; h<K; h++)
         {  
	       bit=( (i+g) >> h) & 1; // get the value of the h-bit of i
	       if(bit)
	        { vecPattern[h]=1; }
	       else
	        { vecPattern[h]=0; }
	     }
	     for(j=0; j<K+1; j++)
	     {
		   if(distances[j]==maxDist)
		    { distances[j]=dist; dummyPointer=Array[j]; Array[j]=dummyArray; dummyArray=dummyPointer; 
		   dummyPointer=ArrayPattern[j]; ArrayPattern[j]=dummyArrayPattern; dummyArrayPattern=dummyPointer; break;}
	     }
	     maxDist=-1;
	     for(j=0; j<K+1; j++)
	     {
		   if(distances[j]>maxDist)
		    { maxDist=distances[j]; }
	     }
	     //mexPrintf("Iteration: %i, New max dist: %1.6f\n",i,maxDist);
	   }
	}
  }

  
  secondpart=clock()-total_start;
  //mexPrintf("Time for part 1: %1.6f, part 2: %1.6f\n",((double)firstpart)/((double)CLOCKS_PER_SEC),((double)secondpart)/((double)CLOCKS_PER_SEC));

  h=0;
  for(i=0;i<K+1;i++)
  {
    // mexPrintf("dist: %f\n",distances[i]);
    distances_out[i]  =  distances[i];
    dummyPointer=Array[i];
	for(j=0;j<M;j++)
	{ Vertices[h]=dummyPointer[j]; h++;} //dummyPointer[j]; h++;}
  } 
  h=0;
  for(i=0;i<K+1;i++)
  {
   dummyPointer=ArrayPattern[i];
	for(j=0;j<K;j++)
	{ VerticesPattern[h]=dummyPointer[j]; h++;} //dummyPointer[j]; h++;}
  }

  delete ArrayGlobal; 
  delete Array;
  delete ArrayPatternGlobal;
  delete ArrayPattern;
  delete PreCompVals;
  delete distances;
  delete indices;
  delete diff;
}
