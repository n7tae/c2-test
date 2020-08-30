/*
 * codec2_fft.h
 *
 *  Created on: 17.09.2016
 *      Author: danilo
 */

#ifndef DRIVERS_FREEDV_CODEC2_FFT_H_
#define DRIVERS_FREEDV_CODEC2_FFT_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

#include "defines.h"

#include "kiss_fftr.h"
#include "kiss_fft.h"

static inline void codec2_fftr(kiss_fftr_state *cfg, float *in, std::complex<float> *out)
{
    kiss_fftr(cfg, in, out);
}

static inline void codec2_fftri(kiss_fftr_state *cfg, std::complex<float> *in, float *out)
{
    kiss_fftri(cfg, in, out);
}

kiss_fft_state *codec2_fft_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem);
kiss_fftr_state *codec2_fftr_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem);
void codec2_fft_free(kiss_fft_state *cfg);
void codec2_fftr_free(kiss_fftr_state *cfg);


static inline void codec2_fft(kiss_fft_state *cfg, std::complex<float> *in, std::complex<float> *out)
{
      kiss_fft(cfg, in, out);
}

void codec2_fft_inplace(kiss_fft_state *cfg, std::complex<float> *inout);

#endif
