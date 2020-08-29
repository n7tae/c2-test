/*---------------------------------------------------------------------------*\

  FILE........: quantise.h
  AUTHOR......: David Rowe
  DATE CREATED: 31/5/92

  Quantisation functions for the sinusoidal coder.

\*---------------------------------------------------------------------------*/

/*
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

#ifndef __QUANTISE__
#define __QUANTISE__

#include "codec2_fft.h"
#include "comp.h"
#include "qbase.h"

class CQuantize : public CQbase {
public:
	void aks_to_M2(kiss_fftr_cfg fftr_fwd_cfg, float ak[], int order, MODEL *model, float E, float *snr, int dump, int sim_pf, int pf, int bass_boost, float beta, float gamma, COMP Aw[]);

	int   encode_Wo(C2CONST *c2const, float Wo, int bits);
	float decode_Wo(C2CONST *c2const, int index, int bits);
	void  encode_lsps_scalar(int indexes[], float lsp[], int order);
	void  decode_lsps_scalar(float lsp[], int indexes[], int order);
	void  encode_lspds_scalar(int indexes[], float lsp[], int order);
	void  decode_lspds_scalar(float lsp[], int indexes[], int order);

	void encode_lsps_vq(int *indexes, float *x, float *xq, int order);
	void decode_lsps_vq(int *indexes, float *xq, int order, int stages);

	int encode_energy(float e, int bits);
	float decode_energy(int index, int bits);

	void pack(unsigned char * bits, unsigned int *nbit, int index, unsigned int index_bits);
	void pack_natural_or_gray(unsigned char * bits, unsigned int *nbit, int index, unsigned int index_bits, unsigned int gray);
	int  unpack(const unsigned char * bits, unsigned int *nbit, unsigned int index_bits);
	int  unpack_natural_or_gray(const unsigned char * bits, unsigned int *nbit, unsigned int index_bits, unsigned int gray);

	int lsp_bits(int i);
	int lspd_bits(int i);
	int lsp_pred_vq_bits(int i);

	void apply_lpc_correction(MODEL *model);
	float speech_to_uq_lsps(float lsp[], float ak[], float Sn[], float w[], int m_pitch, int order);
	int check_lsp_order(float lsp[], int lpc_order);
	void bw_expand_lsps(float lsp[], int order, float min_sep_low, float min_sep_high);

private:
	void compute_weights(const float *x, float *w, int ndim);
	int find_nearest(const float *codebook, int nb_entries, float *x, int ndim);
	void force_min_lsp_dist(float lsp[], int order);
	void lpc_post_filter(kiss_fftr_cfg fftr_fwd_cfg, float Pw[], float ak[], int order, int dump, float beta, float gamma, int bass_boost, float E);
};

#endif
