#ifndef NEWAMPBASE_H
#define NEWQMPBASE_H

#include <complex>
#include "qbase.h"
#include "kiss_fft.h"

class CNewampbase : public CQbase
{
protected:
	void interp_para(float y[], float xp[], float yp[], int np, float x[], int n);
	void resample_const_rate_f(C2CONST *c2const, MODEL *model, float rate_K_vec[], float rate_K_sample_freqs_kHz[], int K);
	void post_filter_newamp1(float vec[], float sample_freq_kHz[], int K, float pf_gain);
	void interp_Wo_v(float Wo_[], int L_[], int voicing_[], float Wo1, float Wo2, int voicing1, int voicing2);
	void determine_phase(C2CONST *c2const, std::complex<float> H[], MODEL *model, int Nfft, kiss_fft_state *fwd_cfg, kiss_fft_state *inv_cfg);
private:
	void mag_to_phase(float phase[], float Gdbfk[], int Nfft, kiss_fft_state *fwd_cfg, kiss_fft_state *inv_cfg);
};

#endif
