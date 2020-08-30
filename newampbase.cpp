#include <assert.h>
#include <math.h>

#include "newampbase.h"

/*---------------------------------------------------------------------------*\

  FUNCTION....: interp_para()
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  General 2nd order parabolic interpolator.  Used splines orginally,
  but this is much simpler and we don't need much accuracy.  Given two
  vectors of points xp and yp, find interpolated values y at points x.

\*---------------------------------------------------------------------------*/

void CNewampbase::interp_para(float y[], float xp[], float yp[], int np, float x[], int n)
{
	assert(np >= 3);

	int k,i;
	float xi, x1, y1, x2, y2, x3, y3, a, b;

	k = 0;
	for (i=0; i<n; i++)
	{
		xi = x[i];

		/* k is index into xp of where we start 3 points used to form parabola */

		while ((xp[k+1] < xi) && (k < (np-3)))
			k++;

		x1 = xp[k];
		y1 = yp[k];
		x2 = xp[k+1];
		y2 = yp[k+1];
		x3 = xp[k+2];
		y3 = yp[k+2];

		//printf("k: %d np: %d i: %d xi: %f x1: %f y1: %f\n", k, np, i, xi, x1, y1);

		a = ((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
		b = ((y3-y2)/(x3-x2)*(x2-x1)+(y2-y1)/(x2-x1)*(x3-x2))/(x3-x1);

		y[i] = a*(xi-x2)*(xi-x2) + b*(xi-x2) + y2;
	}
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: resample_const_rate_f()
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Resample Am from time-varying rate L=floor(pi/Wo) to fixed rate K.

\*---------------------------------------------------------------------------*/

void CNewampbase::resample_const_rate_f(C2CONST *c2const, MODEL *model, float rate_K_vec[], float rate_K_sample_freqs_kHz[], int K)
{
	int m;
	float AmdB[MAX_AMP+1], rate_L_sample_freqs_kHz[MAX_AMP+1], AmdB_peak;

	/* convert rate L=pi/Wo amplitude samples to fixed rate K */

	AmdB_peak = -100.0;
	for(m=1; m<=model->L; m++)
	{
		AmdB[m] = 20.0*log10f(model->A[m]+1E-16);
		if (AmdB[m] > AmdB_peak)
		{
			AmdB_peak = AmdB[m];
		}
		rate_L_sample_freqs_kHz[m] = m*model->Wo*(c2const->Fs/2000.0)/M_PI;
		//printf("m: %d AmdB: %f AmdB_peak: %f  sf: %f\n", m, AmdB[m], AmdB_peak, rate_L_sample_freqs_kHz[m]);
	}

	/* clip between peak and peak -50dB, to reduce dynamic range */

	for(m=1; m<=model->L; m++)
	{
		if (AmdB[m] < (AmdB_peak-50.0))
		{
			AmdB[m] = AmdB_peak-50.0;
		}
	}

	interp_para(rate_K_vec, &rate_L_sample_freqs_kHz[1], &AmdB[1], model->L, rate_K_sample_freqs_kHz, K);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: post_filter
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Post Filter, has a big impact on speech quality after VQ.  When used
  on a mean removed rate K vector, it raises formants, and supresses
  anti-formants.  As it manipulates amplitudes, we normalise energy to
  prevent clipping or large level variations.  pf_gain of 1.2 to 1.5
  (dB) seems to work OK.  Good area for further investigations and
  improvements in speech quality.

\*---------------------------------------------------------------------------*/

void CNewampbase::post_filter_newamp1(float vec[], float sample_freq_kHz[], int K, float pf_gain)
{
	int k;

	/*
	  vec is rate K vector describing spectrum of current frame lets
	  pre-emp before applying PF. 20dB/dec over 300Hz.  Postfilter
	  affects energy of frame so we measure energy before and after
	  and normalise.  Plenty of room for experiment here as well.
	*/

	float pre[K];
	float e_before = 0.0;
	float e_after = 0.0;
	for(k=0; k<K; k++)
	{
		pre[k] = 20.0*log10f(sample_freq_kHz[k]/0.3);
		vec[k] += pre[k];
		e_before += exp10f(vec[k]/10.0);
		vec[k] *= pf_gain;
		e_after += exp10f(vec[k]/10.0);
	}

	float gain = e_after/e_before;
	float gaindB = 10*log10f(gain);

	for(k=0; k<K; k++)
	{
		vec[k] -= gaindB;
		vec[k] -= pre[k];
	}
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: interp_Wo_v
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Decoder side interpolation of Wo and voicing, to go from 25 Hz
  sample rate used over channle to 100Hz internal sample rate of Codec 2.

\*---------------------------------------------------------------------------*/

void CNewampbase::interp_Wo_v(float Wo_[], int L_[], int voicing_[], float Wo1, float Wo2, int voicing1, int voicing2)
{
	int i;
	int M = 4;  /* interpolation rate */

	for(i=0; i<M; i++)
		voicing_[i] = 0;

	if (!voicing1 && !voicing2)
	{
		for(i=0; i<M; i++)
			Wo_[i] = 2.0*M_PI/100.0;
	}

	if (voicing1 && !voicing2)
	{
		Wo_[0] = Wo_[1] = Wo1;
		Wo_[2] = Wo_[3] = 2.0*M_PI/100.0;
		voicing_[0] = voicing_[1] = 1;
	}

	if (!voicing1 && voicing2)
	{
		Wo_[0] = Wo_[1] = 2.0*M_PI/100.0;
		Wo_[2] = Wo_[3] = Wo2;
		voicing_[2] = voicing_[3] = 1;
	}

	if (voicing1 && voicing2)
	{
		float c;
		for(i=0,c=1.0; i<M; i++,c-=1.0/M)
		{
			Wo_[i] = Wo1*c + Wo2*(1.0-c);
			voicing_[i] = 1;
		}
	}

	for(i=0; i<M; i++)
	{
		L_[i] = floorf(M_PI/Wo_[i]);
	}
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: determine_phase
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Given a magnitude spectrum determine a phase spectrum, used for
  phase synthesis with newamp1.

\*---------------------------------------------------------------------------*/

void CNewampbase::determine_phase(C2CONST *c2const, std::complex<float> H[], MODEL *model, int Nfft, kiss_fft_state *fwd_cfg, kiss_fft_state *inv_cfg)
{
	int i,m,b;
	int Ns = Nfft/2+1;
	float Gdbfk[Ns], sample_freqs_kHz[Ns], phase[Ns];
	float AmdB[MAX_AMP+1], rate_L_sample_freqs_kHz[MAX_AMP+1];

	for(m=1; m<=model->L; m++)
	{
		assert(model->A[m] != 0.0);
		AmdB[m] = 20.0*log10f(model->A[m]);
		rate_L_sample_freqs_kHz[m] = (float)m*model->Wo*(c2const->Fs/2000.0)/M_PI;
	}

	for(i=0; i<Ns; i++)
	{
		sample_freqs_kHz[i] = (c2const->Fs/1000.0)*(float)i/Nfft;
	}

	interp_para(Gdbfk, &rate_L_sample_freqs_kHz[1], &AmdB[1], model->L, sample_freqs_kHz, Ns);
	mag_to_phase(phase, Gdbfk, Nfft, fwd_cfg, inv_cfg);

	for(m=1; m<=model->L; m++)
	{
		b = floorf(0.5+m*model->Wo*Nfft/(2.0*M_PI));
		H[m] = std::polar(1.0f, phase[b]);
	}
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: mag_to_phase
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Algorithm for http://www.dsprelated.com/showcode/20.php ported to C.  See
  also Octave function mag_to_phase.m

  Given a magnitude spectrum in dB, returns a minimum-phase phase
  spectra.

\*---------------------------------------------------------------------------*/

void CNewampbase::mag_to_phase(float phase[], /* Nfft/2+1 output phase samples in radians       */
				  float Gdbfk[],              /* Nfft/2+1 postive freq amplitudes samples in dB */
				  int Nfft,
				  kiss_fft_state *fft_fwd_cfg,
				  kiss_fft_state *fft_inv_cfg
				 )
{
	std::complex<float> Sdb[Nfft], c[Nfft], cf[Nfft], Cf[Nfft];
	int  Ns = Nfft/2+1;
	int  i;

	/* install negative frequency components, 1/Nfft takes into
	   account kiss fft lack of scaling on ifft */

	Sdb[0].real(Gdbfk[0]);
	Sdb[0].imag(0);
	for(i=1; i<Ns; i++)
	{
		Sdb[i] = Sdb[Nfft-i] = std::complex<float>(Gdbfk[i], 0);
	}

	/* compute real cepstrum from log magnitude spectrum */

	codec2_fft(fft_inv_cfg, Sdb, c);
	for(i=0; i<Nfft; i++)
	{
		c[i] /= (float)Nfft;
	}

	/* Fold cepstrum to reflect non-min-phase zeros inside unit circle */

	cf[0] = c[0];
	for(i=1; i<Ns-1; i++)
	{
		auto j = Nfft - i;
		cf[i] += c[j];
	}
	cf[Ns-1] = c[Ns-1];
	for(i=Ns; i<Nfft; i++)
	{
		cf[i].real(0);
		cf[i].imag(0);
	}

	/* Cf = dB_magnitude + j * minimum_phase */

	codec2_fft(fft_fwd_cfg, cf, Cf);

	/*  The maths says we are meant to be using log(x), not 20*log10(x),
	    so we need to scale the phase to account for this:
	    log(x) = 20*log10(x)/scale */

	float scale = (20.0/logf(10.0));

	for(i=0; i<Ns; i++)
	{
		phase[i] = Cf[i].imag() / scale;
	}
}
