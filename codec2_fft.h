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

#include "defines.h"
#include "comp.h"


typedef COMP    codec2_fft_cpx;
#include "kiss_fftr.h"

#include "kiss_fft.h"
typedef kiss_fftr_cfg codec2_fftr_cfg;
typedef kiss_fft_cfg codec2_fft_cfg;
typedef kiss_fft_scalar codec2_fft_scalar;



static inline void codec2_fftr(codec2_fftr_cfg cfg, codec2_fft_scalar* in, codec2_fft_cpx* out)
{
    kiss_fftr(cfg, in, (kiss_fft_cpx*)out);
}

static inline void codec2_fftri(codec2_fftr_cfg cfg, codec2_fft_cpx* in, codec2_fft_scalar* out)
{
    kiss_fftri(cfg, (kiss_fft_cpx*)in, out);
}

codec2_fft_cfg codec2_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
codec2_fftr_cfg codec2_fftr_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
void codec2_fft_free(codec2_fft_cfg cfg);
void codec2_fftr_free(codec2_fftr_cfg cfg);


static inline void codec2_fft(codec2_fft_cfg cfg, codec2_fft_cpx* in, codec2_fft_cpx* out)
{
      kiss_fft(cfg, (kiss_fft_cpx*)in, (kiss_fft_cpx*)out);
}

void codec2_fft_inplace(codec2_fft_cfg cfg, codec2_fft_cpx* inout);


#endif
