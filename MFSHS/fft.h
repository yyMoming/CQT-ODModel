#ifndef FFT_H
#define FFT_H
typedef struct
{
	double real;
	double img;
}complex;

typedef struct
{
	float absolute;
	float phase;
}complexWithPhase;


complexWithPhase *fftMain(int fftSize, float *fftBuffer);
complex *fftMain2(int fftSize, double *fftBuffer);
complex *fftMainComplex(int fftSize, complex *fftBuffer);

void fft(int fftSize);
void initW(int fftSize);
void changeFft(int fftSize);
void add(complex, complex, complex *);
void mul(complex, complex, complex *);
void sub(complex, complex, complex *);


void   ifft(int fftSize);
void   divi(complex, complex, complex   *);
complex *ifftMain(int ifftSize, complex *ifftBuffer);
void changeIfft(int fftSize);
complex *fftMain3(int fftSize, float *fftBuffer);
#endif