/*---------------------------------------------------------------------------*\

  FILE........: newamp1.c
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Quantisation functions for the sinusoidal coder, using "newamp1"
  algorithm that resamples variable rate L [Am} to a fixed rate K then
  VQs.

\*---------------------------------------------------------------------------*/

/*
  Copyright David Rowe 2017

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "quantise.h"
#include "mbest.h"
#include "newamp1.h"


/*---------------------------------------------------------------------------*\

  FUNCTION....: ftomel()
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Non linear sampling of frequency axis, reducing the "rate" is a
  first step before VQ

\*---------------------------------------------------------------------------*/

float CNewamp1::ftomel(float fHz)
{
	float mel = floorf(2595.0*log10f(1.0 + fHz/700.0)+0.5);
	return mel;
}

void CNewamp1::mel_sample_freqs_kHz(float rate_K_sample_freqs_kHz[], int K, float mel_start, float mel_end)
{
	float step = (mel_end-mel_start)/(K-1);
	float mel;
	int k;

	mel = mel_start;
	for (k=0; k<K; k++)
	{
		rate_K_sample_freqs_kHz[k] = 0.7*(exp10f(mel/2595.0) - 1.0);
		mel += step;
	}
}



/*---------------------------------------------------------------------------*\

  FUNCTION....: rate_K_mbest_encode
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Two stage rate K newamp1 VQ quantiser using mbest search.

\*---------------------------------------------------------------------------*/

float CNewamp1::rate_K_mbest_encode(int *indexes, float *x, float *xq, int ndim, int mbest_entries)
{
	int i, j, n1, n2;
	const float *codebook1 = newamp1vq_cb[0].cb;
	const float *codebook2 = newamp1vq_cb[1].cb;
	struct MBEST *mbest_stage1, *mbest_stage2;
	float target[ndim];
	float w[ndim];
	int   index[MBEST_STAGES];
	float mse, tmp;

	/* codebook is compiled for a fixed K */

	assert(ndim == newamp1vq_cb[0].k);

	/* equal weights, could be argued mel freq axis gives freq dep weighting */

	for(i=0; i<ndim; i++)
		w[i] = 1.0;

	mbest_stage1 = mbest_create(mbest_entries);
	mbest_stage2 = mbest_create(mbest_entries);
	for(i=0; i<MBEST_STAGES; i++)
		index[i] = 0;

	/* Stage 1 */

	mbest_search(codebook1, x, w, ndim, newamp1vq_cb[0].m, mbest_stage1, index);

	/* Stage 2 */

	for (j=0; j<mbest_entries; j++)
	{
		index[1] = n1 = mbest_stage1->list[j].index[0];
		for(i=0; i<ndim; i++)
			target[i] = x[i] - codebook1[ndim*n1+i];
		mbest_search(codebook2, target, w, ndim, newamp1vq_cb[1].m, mbest_stage2, index);
	}

	n1 = mbest_stage2->list[0].index[1];
	n2 = mbest_stage2->list[0].index[0];
	mse = 0.0;
	for (i=0; i<ndim; i++)
	{
		tmp = codebook1[ndim*n1+i] + codebook2[ndim*n2+i];
		mse += (x[i]-tmp)*(x[i]-tmp);
		xq[i] = tmp;
	}

	mbest_destroy(mbest_stage1);
	mbest_destroy(mbest_stage2);

	indexes[0] = n1;
	indexes[1] = n2;

	return mse;
}




/*---------------------------------------------------------------------------*\

  FUNCTION....: resample_rate_L
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Decoder side conversion of rate K vector back to rate L.

\*---------------------------------------------------------------------------*/

void CNewamp1::resample_rate_L(C2CONST *c2const, MODEL *model, float rate_K_vec[], float rate_K_sample_freqs_kHz[], int K)
{
	float rate_K_vec_term[K+2], rate_K_sample_freqs_kHz_term[K+2];
	float AmdB[MAX_AMP+1], rate_L_sample_freqs_kHz[MAX_AMP+1];
	int m,k;

	/* terminate either end of the rate K vecs with 0dB points */

	rate_K_vec_term[0] = rate_K_vec_term[K+1] = 0.0;
	rate_K_sample_freqs_kHz_term[0] = 0.0;
	rate_K_sample_freqs_kHz_term[K+1] = 4.0;

	for(k=0; k<K; k++)
	{
		rate_K_vec_term[k+1] = rate_K_vec[k];
		rate_K_sample_freqs_kHz_term[k+1] = rate_K_sample_freqs_kHz[k];
		//printf("k: %d f: %f rate_K: %f\n", k, rate_K_sample_freqs_kHz[k], rate_K_vec[k]);
	}

	for(m=1; m<=model->L; m++)
	{
		rate_L_sample_freqs_kHz[m] = m*model->Wo*(c2const->Fs/2000.0)/M_PI;
	}

	interp_para(&AmdB[1], rate_K_sample_freqs_kHz_term, rate_K_vec_term, K+2, &rate_L_sample_freqs_kHz[1], model->L);
	for(m=1; m<=model->L; m++)
	{
		model->A[m] = exp10f(AmdB[m]/20.0);
		// printf("m: %d f: %f AdB: %f A: %f\n", m, rate_L_sample_freqs_kHz[m], AmdB[m], model->A[m]);
	}
}


/* update and optionally run "front eq" equaliser on before VQ */
void CNewamp1::newamp1_eq(float rate_K_vec_no_mean[], float eq[], int K, int eq_en)
{
	static float ideal[] = {8,10,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,-20};
	float gain = 0.02;
	float update;

	for(int k=0; k<K; k++)
	{
		update = rate_K_vec_no_mean[k] - ideal[k];
		eq[k] = (1.0-gain)*eq[k] + gain*update;
		if (eq[k] < 0.0) eq[k] = 0.0;
		if (eq_en)
			rate_K_vec_no_mean[k] -= eq[k];
	}
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: newamp1_model_to_indexes
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  newamp1 encoder for amplitdues {Am}.  Given the rate L model
  parameters, outputs VQ and energy quantiser indexes.

\*---------------------------------------------------------------------------*/

void CNewamp1::newamp1_model_to_indexes(C2CONST *c2const, int indexes[], MODEL *model, float rate_K_vec[], float rate_K_sample_freqs_kHz[], int K, float *mean, float rate_K_vec_no_mean[], float  rate_K_vec_no_mean_[], float *se, float *eq, int eq_en)
{
	int k;

	/* convert variable rate L to fixed rate K */
	resample_const_rate_f(c2const, model, rate_K_vec, rate_K_sample_freqs_kHz, K);

	/* remove mean */
	float sum = 0.0;
	for(k=0; k<K; k++)
		sum += rate_K_vec[k];
	*mean = sum/K;
	for(k=0; k<K; k++)
		rate_K_vec_no_mean[k] = rate_K_vec[k] - *mean;

	/* update and optionally run "front eq" equaliser on before VQ */
	newamp1_eq(rate_K_vec_no_mean, eq, K, eq_en);

	/* two stage VQ */
	rate_K_mbest_encode(indexes, rate_K_vec_no_mean, rate_K_vec_no_mean_, K, NEWAMP1_VQ_MBEST_DEPTH);

	/* running sum of squared error for variance calculation */
	for(k=0; k<K; k++)
		*se += pow(rate_K_vec_no_mean[k]-rate_K_vec_no_mean_[k],2.0);

	/* scalar quantise mean (effectively the frame energy) */
	float w[1] = {1.0};
	float se_mean;
	indexes[2] = quantise(newamp1_energy_cb[0].cb, mean, w, newamp1_energy_cb[0].k, newamp1_energy_cb[0].m, &se_mean);

	/* scalar quantise Wo.  We steal the smallest Wo index to signal
	   an unvoiced frame */
	if (model->voiced)
	{
		int index = encode_log_Wo(c2const, model->Wo, 6);
		if (index == 0)
		{
			index = 1;
		}
		indexes[3] = index;
	}
	else
	{
		indexes[3] = 0;
	}
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: newamp1_interpolate
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

\*---------------------------------------------------------------------------*/

void CNewamp1::newamp1_interpolate(float interpolated_surface_[], float left_vec[], float right_vec[], int K)
{
	int  i, k;
	int  M = 4;
	float c;

	/* (linearly) interpolate 25Hz amplitude vectors back to 100Hz */

	for(i=0,c=1.0; i<M; i++,c-=1.0/M)
	{
		for(k=0; k<K; k++)
		{
			interpolated_surface_[i*K+k] = left_vec[k]*c + right_vec[k]*(1.0-c);
		}
	}
}


/*---------------------------------------------------------------------------*\

  FUNCTION....: newamp1_indexes_to_rate_K_vec
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  newamp1 decoder for amplitudes {Am}.  Given the rate K VQ and energy
  indexes, outputs rate K vector.

\*---------------------------------------------------------------------------*/

void CNewamp1::newamp1_indexes_to_rate_K_vec(float rate_K_vec_[], float rate_K_vec_no_mean_[], float rate_K_sample_freqs_kHz[], int K, float *mean_, int indexes[], float user_rate_K_vec_no_mean_[], int post_filter_en)
{
	int   k;
	const float *codebook1 = newamp1vq_cb[0].cb;
	const float *codebook2 = newamp1vq_cb[1].cb;
	int n1 = indexes[0];
	int n2 = indexes[1];

	if (user_rate_K_vec_no_mean_ == NULL)
	{
		/* normal operation */
		for(k=0; k<K; k++)
		{
			rate_K_vec_no_mean_[k] = codebook1[K*n1+k] + codebook2[K*n2+k];
		}
	}
	else
	{
		/* for development we can optionally inject the quantised rate K vector here */
		for(k=0; k<K; k++)
			rate_K_vec_no_mean_[k] = user_rate_K_vec_no_mean_[k];
	}

	if (post_filter_en)
		post_filter_newamp1(rate_K_vec_no_mean_, rate_K_sample_freqs_kHz, K, 1.5);

	*mean_ = newamp1_energy_cb[0].cb[indexes[2]];

	for(k=0; k<K; k++)
	{
		rate_K_vec_[k] = rate_K_vec_no_mean_[k] + *mean_;
	}
}


/*---------------------------------------------------------------------------*\

  FUNCTION....: newamp1_indexes_to_model
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  newamp1 decoder.

\*---------------------------------------------------------------------------*/

void CNewamp1::newamp1_indexes_to_model(C2CONST *c2const, MODEL model_[], std::complex<float> H[], float *interpolated_surface_, float  prev_rate_K_vec_[], float *Wo_left, int *voicing_left, float  rate_K_sample_freqs_kHz[], int K, kiss_fft_state *fwd_cfg, kiss_fft_state *inv_cfg, int indexes[], float  user_rate_K_vec_no_mean_[], int post_filter_en)
{
	float rate_K_vec_[K], rate_K_vec_no_mean_[K], mean_, Wo_right;
	int   voicing_right, k;
	int   M = 4;

	/* extract latest rate K vector */

	newamp1_indexes_to_rate_K_vec(rate_K_vec_, rate_K_vec_no_mean_, rate_K_sample_freqs_kHz, K, &mean_, indexes, user_rate_K_vec_no_mean_, post_filter_en);

	/* decode latest Wo and voicing */

	if (indexes[3])
	{
		Wo_right = decode_log_Wo(c2const, indexes[3], 6);
		voicing_right = 1;
	}
	else
	{
		Wo_right  = 2.0*M_PI/100.0;
		voicing_right = 0;
	}

	/* interpolate 25Hz rate K vec back to 100Hz */

	float *left_vec = prev_rate_K_vec_;
	float *right_vec = rate_K_vec_;
	newamp1_interpolate(interpolated_surface_, left_vec, right_vec, K);

	/* interpolate 25Hz v and Wo back to 100Hz */

	float aWo_[M];
	int avoicing_[M], aL_[M], i;

	interp_Wo_v(aWo_, aL_, avoicing_, *Wo_left, Wo_right, *voicing_left, voicing_right);

	/* back to rate L amplitudes, Synthesise phase for each frame */

	for(i=0; i<M; i++)
	{
		model_[i].Wo = aWo_[i];
		model_[i].L  = aL_[i];
		model_[i].voiced = avoicing_[i];

		resample_rate_L(c2const, &model_[i], &interpolated_surface_[K*i], rate_K_sample_freqs_kHz, K);
		determine_phase(c2const, &H[(MAX_AMP+1)*i], &model_[i], NEWAMP1_PHASE_NFFT, fwd_cfg, inv_cfg);
	}

	/* update memories for next time */

	for(k=0; k<K; k++)
	{
		prev_rate_K_vec_[k] = rate_K_vec_[k];
	}
	*Wo_left = Wo_right;
	*voicing_left = voicing_right;

}
