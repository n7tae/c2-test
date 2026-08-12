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
	FILE          *fin;
	FILE          *fout;

	if (argc < 3)
	{
		printf("usage: c2enc InputRawspeechFile OutputBitFile\n");
		exit(1);
	}

	if (strcmp(argv[1], "-")  == 0) fin = stdin;
	else if ( (fin = fopen(argv[1],"rb")) == NULL )
	{
		fprintf(stderr, "Error opening input speech file: %s: %s.\n",
				argv[1], strerror(errno));
		exit(1);
	}

	if (strcmp(argv[2], "-") == 0) fout = stdout;
	else if ( (fout = fopen(argv[2],"wb")) == NULL )
	{
		fprintf(stderr, "Error opening output compressed bit file: %s: %s.\n",
				argv[2], strerror(errno));
		exit(1);
	}

	CCodec2 cc2;
	auto nsam = cc2.samples_per_frame();
	auto nbit = cc2.bits_per_frame();
	short buf[nsam];
	auto nbyte = (nbit + 7) / 8;

	unsigned char bits[nbyte];

	//fprintf(stderr,"gray: %d softdec: %d\n", gray, softdec);

	while(fread(buf, sizeof(short), nsam, fin) == (size_t)nsam)
	{

		cc2.encode(bits, buf);

		fwrite(bits, sizeof(char), nbyte, fout);

		// if this is in a pipeline, we probably don't want the usual
		// buffering to occur

		if (fout == stdout) fflush(stdout);
		if (fin == stdin) fflush(stdin);
	}

	fclose(fin);
	fclose(fout);

	return 0;
}
