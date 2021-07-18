#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include <stdio.h>


complex *ifftTransform;
complex *fftTransform, *W;
float PI;


complex *fftMainComplex(int fftSize, complex *fftBuffer)
{

    fftTransform = (complex *)malloc(sizeof(complex)* fftSize);
    PI = atan(1.0) * 4;

//    complexWithPhase *fftTransformWithPhase;
//    fftTransformWithPhase = (complexWithPhase *)malloc(sizeof(complexWithPhase)* fftSize);

    complex *fftBufferMiddleVar;
    fftBufferMiddleVar = fftBuffer;
    for (int i = 0; i<fftSize; i++, fftBufferMiddleVar++)
    {
        fftTransform[i].real = (*fftBufferMiddleVar).real;
        fftTransform[i].img = (*fftBufferMiddleVar).img;
    }

    initW(fftSize);
    fft(fftSize);


    free(W);
    return fftTransform;
}

complexWithPhase *fftMain(int fftSize, float *fftBuffer)
{	
	fftTransform = (complex *)malloc(sizeof(complex)* fftSize);
	PI = atan(1.0) * 4;

	complexWithPhase *fftTransformWithPhase;
	fftTransformWithPhase = (complexWithPhase *)malloc(sizeof(complexWithPhase)* fftSize);
	float *fftBufferMiddleVar;
	fftBufferMiddleVar = fftBuffer;
	for (int i = 0; i<fftSize; i++, fftBufferMiddleVar++)
	{
		fftTransform[i].real = *fftBufferMiddleVar;
		fftTransform[i].img = 0;
	}

	initW(fftSize);
	fft(fftSize);

	for (int i = 0; i<fftSize; i++)
	{
		fftTransformWithPhase[i].absolute = sqrt(pow(fftTransform[i].real, 2) + pow(fftTransform[i].img, 2));
	}

	for (int i = 0; i<fftSize; i++)
	{
		fftTransformWithPhase[i].phase = atan2(fftTransform[i].img, fftTransform[i].real);
	}
	free(fftTransform);
	free(W);

	
	return fftTransformWithPhase;
}

complex *fftMain2(int fftSize, double *fftBuffer)
{

	fftTransform = (complex *)malloc(sizeof(complex)* fftSize);
	PI = atan(1.0) * 4;

	complexWithPhase *fftTransformWithPhase;
	fftTransformWithPhase = (complexWithPhase *)malloc(sizeof(complexWithPhase)* fftSize);

	double *fftBufferMiddleVar;
	fftBufferMiddleVar = fftBuffer;
	for (int i = 0; i<fftSize; i++, fftBufferMiddleVar++)
	{
		fftTransform[i].real = *fftBufferMiddleVar;
		fftTransform[i].img = 0;
	}

	initW(fftSize);
	fft(fftSize);


	free(W);
	return fftTransform;
}


complex *fftMain3(int fftSize, float *fftBuffer)
{

    fftTransform = (complex *)malloc(sizeof(complex)* fftSize);
    PI = atan(1.0) * 4;

    float *fftBufferMiddleVar;
    fftBufferMiddleVar = fftBuffer;
    for (int i = 0; i<fftSize; i++, fftBufferMiddleVar++)
    {
        fftTransform[i].real = *fftBufferMiddleVar;
        fftTransform[i].img = 0;
    }

    initW(fftSize);
    fft(fftSize);


    free(W);
    return fftTransform;
}
complex *ifftMain(int ifftSize, complex *ifftBuffer)
{
	ifftTransform = (complex *)malloc(sizeof(complex)* ifftSize);
	PI = 3.141592653589793;//atan(1.0)*4;

	for (int i = 0; i<ifftSize; i++)
	{
		ifftTransform[i].real = ifftBuffer[i].real;
		ifftTransform[i].img = ifftBuffer[i].img;
	}

	initW(ifftSize);
	ifft(ifftSize);

	free(W);
	W = NULL;

	return ifftTransform;
}

void ifft(int ifftSize)
{
	int   i = 0, j = 0, k = 0, l = ifftSize;
	complex   up, down;
	for (i = 0; i< log((double)ifftSize) / log((double)2); i++)
	{
		l /= 2;
		for (j = 0; j<ifftSize; j += 2 * l)
		{
			for (k = 0; k<l; k++)
			{
				add(ifftTransform[j + k], ifftTransform[j + k + l], &up);
				up.real /= 2; up.img /= 2;
				sub(ifftTransform[j + k], ifftTransform[j + k + l], &down);
				down.real /= 2; down.img /= 2;
				divi(down, W[ifftSize*k / 2 / l], &down);
				ifftTransform[j + k] = up;
				ifftTransform[j + k + l] = down;
			}
		}
	}
	changeIfft(ifftSize);
}

void fft(int fftSize)
{
	int i = 0, j = 0, k = 0, l = 0;
	complex up, down, product;
	changeFft(fftSize);
	for (i = 0; i< log((double)fftSize) / log((double)2); i++){
		l = 1 << i;
		for (j = 0; j<fftSize; j += 2 * l){
			for (k = 0; k<l; k++){
				mul(fftTransform[j + k + l], W[fftSize*k / 2 / l], &product);
				add(fftTransform[j + k], product, &up);
				sub(fftTransform[j + k], product, &down);
				fftTransform[j + k] = up;
				fftTransform[j + k + l] = down;
			}
		}
	}
}


void initW(int fftSize)
{
	int i;
	W = (complex *)malloc(sizeof(complex)* fftSize);
	for (i = 0; i<fftSize; i++){
		W[i].real = cos(2 * PI / fftSize*i);
		W[i].img = -1 * sin(2 * PI / fftSize*i);
	}
}


void changeFft(int fftSize)
{
	complex temp;
	unsigned short i = 0, j = 0, k = 0;
	double t;
	for (i = 0; i<fftSize; i++){
		k = i; j = 0;
		t = log((double)fftSize) / log((double)2);
		while ((t--)>0){
			j = j << 1;
			j |= (k & 1);
			k = k >> 1;
		}
		if (j>i){
			temp = fftTransform[i];
			fftTransform[i] = fftTransform[j];
			fftTransform[j] = temp;
		}
	}
}

void   changeIfft(int ifftSize)
{
	complex   temp;
	unsigned   short   i = 0, j = 0, k = 0;
	double   t;
	for (i = 0; i<ifftSize; i++)
	{
		k = i; j = 0;
		t = (log((double)ifftSize)) / log((double)2);
		while ((t--)>0)
		{
			j = j << 1;
			j |= (k & 1);
			k = k >> 1;
		}
		if (j>i)
		{
			temp = ifftTransform[i];
			ifftTransform[i] = ifftTransform[j];
			ifftTransform[j] = temp;
		}
	}
}


void add(complex a, complex b, complex *c)        
{
	c->real = a.real + b.real;
	c->img = a.img + b.img;
}

void mul(complex a, complex b, complex *c)
{
	c->real = a.real*b.real - a.img*b.img;
	c->img = a.real*b.img + a.img*b.real;
}
void sub(complex a, complex b, complex *c)
{
	c->real = a.real - b.real;
	c->img = a.img - b.img;
}

void  divi(complex   a, complex   b, complex   *c)
{
	c->real = (a.real*b.real + a.img*b.img) / (
		b.real*b.real + b.img*b.img);
	c->img = (a.img*b.real - a.real*b.img) / (b.real*b.real + b.img*b.img);
}
