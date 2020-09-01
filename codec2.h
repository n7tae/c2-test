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

#include "codec2_internal.h"
#include "defines.h"
#include "kiss_fft.h"
#include "version.h"
#include "nlp.h"
#include "quantise.h"
#include "newamp1.h"
#include "newamp2.h"

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

class CCodec2
{
public:
	bool codec2_create(int mode);
	void codec2_destroy();
	void codec2_encode(unsigned char * bits, short speech_in[]);
	void codec2_decode(short speech_out[], const unsigned char *bits);
	int  codec2_samples_per_frame();
	int  codec2_bits_per_frame();
	void codec2_open_mlfeat(char *feat_filename, char *model_filename);
	void codec2_load_codebook(int num, char *filename);
	void codec2_set_natural_or_gray(int gray);
	void codec2_700c_eq(int en);
	float codec2_get_var();
	float *codec2_enable_user_ratek(int *K);
	void codec2_700c_post_filter(int en);
	void codec2_set_softdec(float *softdec);
	float codec2_get_energy(const unsigned char *bits);
	void codec2_decode_ber(short speech_out[], const unsigned char *bits, float ber_est);

private:
	void codec2_set_lpc_post_filter(int enable, int bass_boost, float beta, float gamma);
	int  codec2_get_spare_bit_index();
	int  codec2_rebuild_spare_bit(char unpacked_bits[]);

	// merged from other files
	void sample_phase(MODEL *model, std::complex<float> filter_phase[], std::complex<float> A[]);
	void phase_synth_zero_order(int n_samp, MODEL *model, float *ex_phase, std::complex<float> filter_phase[]);
	void postfilter(MODEL *model, float *bg_est);

	C2CONST c2const_create(int Fs, float framelength_ms);

	void make_analysis_window(C2CONST *c2const, FFT_STATE *fft_fwd_cfg, float w[], float W[]);
	void dft_speech(C2CONST *c2const, FFT_STATE &fft_fwd_cfg, std::complex<float> Sw[], float Sn[], float w[]);
	void two_stage_pitch_refinement(C2CONST *c2const, MODEL *model, std::complex<float> Sw[]);
	void estimate_amplitudes(MODEL *model, std::complex<float> Sw[], float W[], int est_phase);
	float est_voicing_mbe(C2CONST *c2const, MODEL *model, std::complex<float> Sw[], float W[]);
	void make_synthesis_window(C2CONST *c2const, float Pn[]);
	void synthesise(int n_samp, FFTR_STATE *fftr_inv_cfg, float Sn_[], MODEL *model, float Pn[], int shift);
	int codec2_rand(void);
	void hs_pitch_refinement(MODEL *model, std::complex<float> Sw[], float pmin, float pmax, float pstep);

	void interp_Wo(MODEL *interp, MODEL *prev, MODEL *next, float Wo_min);
	void interp_Wo2(MODEL *interp, MODEL *prev, MODEL *next, float weight, float Wo_min);
	float interp_energy(float prev, float next);
	float interp_energy2(float prev, float next, float weight);
	void interpolate_lsp_ver2(float interp[], float prev[],  float next[], float weight, int order);

	void analyse_one_frame(MODEL *model, short speech[]);
	void synthesise_one_frame(short speech[], MODEL *model, std::complex<float> Aw[], float gain);
	void codec2_encode_3200(unsigned char * bits, short speech[]);
	void codec2_decode_3200(short speech[], const unsigned char * bits);
	void codec2_encode_2400(unsigned char * bits, short speech[]);
	void codec2_decode_2400(short speech[], const unsigned char * bits);
	void codec2_encode_1600(unsigned char * bits, short speech[]);
	void codec2_decode_1600(short speech[], const unsigned char * bits);
	void codec2_encode_1400(unsigned char * bits, short speech[]);
	void codec2_decode_1400(short speech[], const unsigned char * bits);
	void codec2_encode_1300(unsigned char * bits, short speech[]);
	void codec2_decode_1300(short speech[], const unsigned char * bits, float ber_est);
	void codec2_encode_1200(unsigned char * bits, short speech[]);
	void codec2_decode_1200(short speech[], const unsigned char * bits);
	void codec2_encode_700c(unsigned char * bits, short speech[]);
	void codec2_decode_700c(short speech[], const unsigned char * bits);
	void codec2_encode_450(unsigned char * bits, short speech[]);
	void codec2_decode_450(short speech[], const unsigned char * bits);
	void codec2_decode_450pwb(short speech[], const unsigned char * bits);
	void ear_protection(float in_out[], int n);
	float codec2_energy_700c(const unsigned char * bits);
	float codec2_energy_450(const unsigned char * bits);

	void (CCodec2::*encode)(unsigned char * bits, short speech[]);
	void (CCodec2::*decode)(short speech[], const unsigned char * bits);
	void (CCodec2::*decode_ber)(short speech[], const unsigned char * bits, float ber_est);
	Cnlp nlp;
	CQuantize qt;
	CNewamp1 na1;
	CNewamp2 na2;
	CODEC2 c2;
};

#endif
