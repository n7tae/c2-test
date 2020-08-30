/*
Copyright (c) 2003-2004, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "defines.h"
#include "kiss_fftr.h"
#include "assert.h"

kiss_fftr_state *kiss_fftr_alloc(int nfft, int inverse_fft, void * mem, size_t *lenmem)
{
	int i;
	kiss_fftr_state *st = NULL;
	size_t subsize, memneeded;

	if (nfft & 1)
	{
		fprintf(stderr,"Real FFT optimization must be even.\n");
		return NULL;
	}
	nfft >>= 1;

	kiss_fft_alloc (nfft, inverse_fft, NULL, &subsize);
	memneeded = sizeof(struct kiss_fftr_state) + subsize + sizeof(kiss_fft_cpx) * ( nfft * 3 / 2);

	if (lenmem == NULL)
	{
		st = (kiss_fftr_state *)malloc(memneeded);
	}
	else
	{
		if (*lenmem >= memneeded)
			st = (kiss_fftr_state *)mem;
		*lenmem = memneeded;
	}
	if (!st)
		return NULL;

	st->substate = (kiss_fft_state *)(st + 1); /*just beyond kiss_fftr_state struct */
	st->tmpbuf = (kiss_fft_cpx *) (((char *) st->substate) + subsize);
	st->super_twiddles = st->tmpbuf + nfft;
	kiss_fft_alloc(nfft, inverse_fft, st->substate, &subsize);

	for (i = 0; i < nfft/2; ++i)
	{
		float phase =
			-3.14159265358979323846264338327 * ((float) (i+1) / nfft + .5);
		if (inverse_fft)
			phase *= -1;
		st->super_twiddles[i] = std::polar(1.0f, phase);
	}
	return st;
}

void kiss_fftr(kiss_fftr_state *st, const float *timedata, kiss_fft_cpx *freqdata)
{
	/* input buffer timedata is stored row-wise */
	int k,ncfft;
	kiss_fft_cpx fpnk,fpk,f1k,f2k,tw,tdc;

	assert(st->substate->inverse==0);

	ncfft = st->substate->nfft;

	/*perform the parallel fft of two real signals packed in real,imag*/
	kiss_fft( st->substate, (const kiss_fft_cpx*)timedata, st->tmpbuf );
	/* The real part of the DC element of the frequency spectrum in st->tmpbuf
	 * contains the sum of the even-numbered elements of the input time sequence
	 * The imag part is the sum of the odd-numbered elements
	 *
	 * The sum of tdc.r and tdc.i is the sum of the input time sequence.
	 *      yielding DC of input time sequence
	 * The difference of tdc.r - tdc.i is the sum of the input (dot product) [1,-1,1,-1...
	 *      yielding Nyquist bin of input time sequence
	 */

	tdc = st->tmpbuf[0];
	freqdata[0].real(tdc.real() + tdc.imag());
	freqdata[ncfft].real(tdc.real() - tdc.imag());
	freqdata[ncfft].imag(0.f);
	freqdata[0].imag(0.f);

	for ( k=1; k <= ncfft/2 ; ++k )
	{
		fpk = st->tmpbuf[k];
		fpnk = std::conj(st->tmpbuf[ncfft-k]);

		f1k = fpk + fpnk;
		f2k = fpk - fpnk;
		tw = f2k * st->super_twiddles[k-1];

		freqdata[k] = 0.5f * (f1k + tw);
		freqdata[ncfft-k].real(0.5f * (f1k.real() - tw.real()));
		freqdata[ncfft-k].imag(0.5f * (tw.imag() - f1k.imag()));
	}
}

void kiss_fftri(kiss_fftr_state *st, const kiss_fft_cpx *freqdata, float *timedata)
{
	/* input buffer timedata is stored row-wise */
	int k, ncfft;

	assert(st->substate->inverse == 1);

	ncfft = st->substate->nfft;

	st->tmpbuf[0].real(freqdata[0].real() + freqdata[ncfft].real());
	st->tmpbuf[0].imag(freqdata[0].real() - freqdata[ncfft].real());

	for (k = 1; k <= ncfft / 2; ++k)
	{
		kiss_fft_cpx fk, fnkc, fek, fok, tmp;
		fk = freqdata[k];
		fnkc = std::conj(freqdata[ncfft - k]);

		fek = fk + fnkc;
		tmp = fk - fnkc;
		fok = tmp * st->super_twiddles[k-1];
		st->tmpbuf[k] = fek + fok;
		st->tmpbuf[ncfft - k] = std::conj(fek - fok);
	}
	kiss_fft (st->substate, st->tmpbuf, (kiss_fft_cpx *) timedata);
}
