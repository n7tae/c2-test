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


#include "kiss_fftr.h"
#include "kiss_fft.h"

static inline void codec2_fftr(kiss_fftr_cfg cfg, float* in, COMP* out)
{
    kiss_fftr(cfg, in, (kiss_fft_cpx*)out);
}

static inline void codec2_fftri(kiss_fftr_cfg cfg, COMP* in, float* out)
{
    kiss_fftri(cfg, (kiss_fft_cpx*)in, out);
}

kiss_fft_cfg codec2_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
kiss_fftr_cfg codec2_fftr_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
void codec2_fft_free(kiss_fft_cfg cfg);
void codec2_fftr_free(kiss_fftr_cfg cfg);


static inline void codec2_fft(kiss_fft_cfg cfg, COMP* in, COMP* out)
{
      kiss_fft(cfg, (kiss_fft_cpx*)in, (kiss_fft_cpx*)out);
}

void codec2_fft_inplace(kiss_fft_cfg cfg, COMP* inout);


#endif
