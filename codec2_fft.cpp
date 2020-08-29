/*
 * codec2_fft.c
 *
 *  Created on: 24.09.2016
 *      Author: danilo
 */

#include "codec2_fft.h"

#include "debug_alloc.h"

#include "_kiss_fft_guts.h"

void codec2_fft_free(kiss_fft_cfg cfg)
{
	free(cfg);
}

kiss_fft_cfg codec2_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem)
{
	kiss_fft_cfg retval;
	retval = kiss_fft_alloc(nfft, inverse_fft, mem, lenmem);
	return retval;
}

kiss_fftr_cfg codec2_fftr_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem)
{
	kiss_fftr_cfg retval;
	retval = kiss_fftr_alloc(nfft, inverse_fft, mem, lenmem);
	return retval;
}
void codec2_fftr_free(kiss_fftr_cfg cfg)
{
	free(cfg);
}

// there is a little overhead for inplace kiss_fft but this is
// on the powerful platforms like the Raspberry or even x86 PC based ones
// not noticeable
// the reduced usage of RAM and increased performance on STM32 platforms
// should be worth it.
void codec2_fft_inplace(kiss_fft_cfg cfg, COMP* inout)
{
	kiss_fft_cpx in[512];
	// decide whether to use the local stack based buffer for in
	// or to allow kiss_fft to allocate RAM
	// second part is just to play safe since first method
	// is much faster and uses less RAM
	if (cfg->nfft <= 512)
	{
		memcpy(in,inout,cfg->nfft*sizeof(kiss_fft_cpx));
		kiss_fft(cfg, in, (kiss_fft_cpx*)inout);
	}
	else
	{
		kiss_fft(cfg, (kiss_fft_cpx*)inout, (kiss_fft_cpx*)inout);
	}
}
