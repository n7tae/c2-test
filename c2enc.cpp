/*---------------------------------------------------------------------------*\

  FILE........: c2enc.c
  AUTHOR......: David Rowe
  DATE CREATED: 23/8/2010

  Encodes a file of raw speech samples using codec2 and outputs a file
  of bits.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "codec2.h"

int main(int argc, char *argv[])
{
	int            mode;
	FILE          *fin;
	FILE          *fout;
	int            bit, byte,i;
	int            eq = 0;

	CCodec2 cc2;

	if (argc < 4)
	{
		printf("usage: c2enc 3200|1600 InputRawspeechFile OutputBitFile\n");
		exit(1);
	}

	if (strcmp(argv[1],"3200") == 0)
		mode = CODEC2_MODE_3200;
	else if (strcmp(argv[1],"1600") == 0)
		mode = CODEC2_MODE_1600;
	else
	{
		fprintf(stderr, "Error in mode: %s.  Must be 3200 or 1600\n", argv[1]);
		exit(1);
	}

	if (strcmp(argv[2], "-")  == 0) fin = stdin;
	else if ( (fin = fopen(argv[2],"rb")) == NULL )
	{
		fprintf(stderr, "Error opening input speech file: %s: %s.\n",
				argv[2], strerror(errno));
		exit(1);
	}

	if (strcmp(argv[3], "-") == 0) fout = stdout;
	else if ( (fout = fopen(argv[3],"wb")) == NULL )
	{
		fprintf(stderr, "Error opening output compressed bit file: %s: %s.\n",
				argv[3], strerror(errno));
		exit(1);
	}

	if (cc2.codec2_create(mode))
		return 1;
	auto nsam = cc2.codec2_samples_per_frame();
	auto nbit = cc2.codec2_bits_per_frame();
	short buf[nsam];
	auto nbyte = (nbit + 7) / 8;

	unsigned char bits[nbyte];

	int gray = 1;
	cc2.codec2_set_natural_or_gray(gray);

	//fprintf(stderr,"gray: %d softdec: %d\n", gray, softdec);

	while(fread(buf, sizeof(short), nsam, fin) == (size_t)nsam)
	{

		cc2.codec2_encode(bits, buf);

		fwrite(bits, sizeof(char), nbyte, fout);

		// if this is in a pipeline, we probably don't want the usual
		// buffering to occur

		if (fout == stdout) fflush(stdout);
		if (fin == stdin) fflush(stdin);
	}

	cc2.codec2_destroy();

	fclose(fin);
	fclose(fout);

	return 0;
}
