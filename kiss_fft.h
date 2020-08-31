#ifndef KISS_FFT_H
#define KISS_FFT_H

#include <complex>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* for real ffts, we need an even size */
#define kiss_fftr_next_fast_size_real(n) (kiss_fft_next_fast_size( ((n)+1) >> 1) << 1 )

class CKissFFT
{
public:
	kiss_fft_state *fft_alloc(int nfft, int inverse_fft, void *mem,size_t *lenmem);
	void fft(kiss_fft_state *cfg, const std::complex<float> *fin, std::complex<float> *fout);
	void fft_stride(kiss_fft_state *cfg, const std::complex<float> *fin, std::complex<float> *fout, int fin_stride);
	int fft_next_fast_size(int n);
	kiss_fftr_state *fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
	void fftr(kiss_fftr_state *cfg,const float *timedata,std::complex<float> *freqdata);
	void fftri(kiss_fftr_state *cfg,const std::complex<float> *freqdata,float *timedata);
private:
	void kf_bfly2(std::complex<float> *Fout, const size_t fstride, kiss_fft_state *st, int m);
	void kf_bfly3(std::complex<float> *Fout, const size_t fstride, kiss_fft_state *st, int m);
	void kf_bfly4(std::complex<float> *Fout, const size_t fstride, kiss_fft_state *st, int m);
	void kf_bfly5(std::complex<float> *Fout, const size_t fstride, kiss_fft_state *st, int m);
	void kf_bfly_generic(std::complex<float> *Fout, const size_t fstride, kiss_fft_state *st, int m, int p);
	void kf_work(std::complex<float> *Fout, const std::complex<float> *f, const size_t fstride, int in_stride, int *factors, kiss_fft_state *st);
	void kf_factor(int n,int * facbuf);
};
#endif
