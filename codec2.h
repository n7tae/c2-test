/*---------------------------------------------------------------------------*\

  FILE........: codec2.h
  AUTHOR......: David Rowe
  DATE CREATED: 21 August 2010

  Codec 2 fully quantised encoder and decoder functions.  If you want use
  Codec 2, these are the functions you need to call.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#ifndef __CODEC2__
#define  __CODEC2__

#include <complex>

#include "defines.h"
#include "kiss_fft.h"
#include "version.h"

#define CODEC2_MODE_3200 	0
#define CODEC2_MODE_2400 	1
#define CODEC2_MODE_1600 	2
#define CODEC2_MODE_1400 	3
#define CODEC2_MODE_1300 	4
#define CODEC2_MODE_1200 	5
#define CODEC2_MODE_700C 	8
#define CODEC2_MODE_450 	10
#define CODEC2_MODE_450PWB 	11

#ifndef CODEC2_MODE_EN_DEFAULT
#define CODEC2_MODE_EN_DEFAULT 1
#endif

#define CODEC2_RAND_MAX 32767

// by default we enable all modes
// disable during compile time with -DCODEC2_MODE_1600_EN=0
// all but CODEC2 1600 are enabled then

//or the other way round
// -DCODEC2_MODE_EN_DEFAULT=0 -DCODEC2_MODE_1600_EN=1
// only CODEC2 Mode 1600

#if !defined(CODEC2_MODE_3200_EN)
        #define CODEC2_MODE_3200_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_2400_EN)
        #define CODEC2_MODE_2400_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1600_EN)
        #define CODEC2_MODE_1600_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1400_EN)
        #define CODEC2_MODE_1400_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1300_EN)
        #define CODEC2_MODE_1300_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1200_EN)
        #define CODEC2_MODE_1200_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_700C_EN)
        #define CODEC2_MODE_700C_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_450_EN)
        #define CODEC2_MODE_450_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_450PWB_EN)
        #define CODEC2_MODE_450PWB_EN CODEC2_MODE_EN_DEFAULT
#endif

#define CODEC2_MODE_ACTIVE(mode_name, var)  ((mode_name##_EN) == 0 ? 0: (var) == mode_name)

struct CODEC2;

struct CODEC2 *codec2_create(int mode);
void codec2_destroy(struct CODEC2 *codec2_state);
void codec2_encode(struct CODEC2 *codec2_state, unsigned char * bits, short speech_in[]);
void codec2_decode(struct CODEC2 *codec2_state, short speech_out[], const unsigned char *bits);
void codec2_decode_ber(struct CODEC2 *codec2_state, short speech_out[], const unsigned char *bits, float ber_est);
int  codec2_samples_per_frame(struct CODEC2 *codec2_state);
int  codec2_bits_per_frame(struct CODEC2 *codec2_state);

void codec2_set_lpc_post_filter(struct CODEC2 *codec2_state, int enable, int bass_boost, float beta, float gamma);
int  codec2_get_spare_bit_index(struct CODEC2 *codec2_state);
int  codec2_rebuild_spare_bit(struct CODEC2 *codec2_state, char unpacked_bits[]);
void codec2_set_natural_or_gray(struct CODEC2 *codec2_state, int gray);
void codec2_set_softdec(struct CODEC2 *c2, float *softdec);
float codec2_get_energy(struct CODEC2 *codec2_state, const unsigned char *bits);

// support for ML and VQ experiments
void codec2_open_mlfeat(struct CODEC2 *codec2_state, char *feat_filename, char *model_filename);
void codec2_load_codebook(struct CODEC2 *codec2_state, int num, char *filename);
float codec2_get_var(struct CODEC2 *codec2_state);
float *codec2_enable_user_ratek(struct CODEC2 *codec2_state, int *K);

// 700C post filter and equaliser
void codec2_700c_post_filter(struct CODEC2 *codec2_state, int en);
void codec2_700c_eq(struct CODEC2 *codec2_state, int en);

// merged from other files
void sample_phase(MODEL *model, std::complex<float> filter_phase[], std::complex<float> A[]);
void phase_synth_zero_order(int n_samp, MODEL *model, float *ex_phase, std::complex<float> filter_phase[]);
void postfilter(MODEL *model, float *bg_est);

C2CONST c2const_create(int Fs, float framelength_ms);

void make_analysis_window(C2CONST *c2const, kiss_fft_state *fft_fwd_cfg, float w[], float W[]);
void dft_speech(C2CONST *c2const, kiss_fft_state *fft_fwd_cfg, std::complex<float> Sw[], float Sn[], float w[]);
void two_stage_pitch_refinement(C2CONST *c2const, MODEL *model, std::complex<float> Sw[]);
void estimate_amplitudes(MODEL *model, std::complex<float> Sw[], float W[], int est_phase);
float est_voicing_mbe(C2CONST *c2const, MODEL *model, std::complex<float> Sw[], float W[]);
void make_synthesis_window(C2CONST *c2const, float Pn[]);
void synthesise(int n_samp, kiss_fftr_state *fftr_inv_cfg, float Sn_[], MODEL *model, float Pn[], int shift);
int codec2_rand(void);
void hs_pitch_refinement(MODEL *model, std::complex<float> Sw[], float pmin, float pmax, float pstep);

void interp_Wo(MODEL *interp, MODEL *prev, MODEL *next, float Wo_min);
void interp_Wo2(MODEL *interp, MODEL *prev, MODEL *next, float weight, float Wo_min);
float interp_energy(float prev, float next);
float interp_energy2(float prev, float next, float weight);
void interpolate_lsp_ver2(float interp[], float prev[],  float next[], float weight, int order);

#endif
