# c2-test

## Description

This is a C++ implementation of David Row's Codec2. The [original source repo](https://github.com/drowe67/codec2) contains the Codec2 vocoder sources and additional code to buid a number of executables. This repo only contains the vocoder sources and source for two test programs, an encoder and a decoder, and this code can only handle 3200 and 1600 modes.

Many changes were performed to translate David's C code to be supported in a C++11 standard. All operational code has been encapsulated into a few classes. The classes selected are based on C++ scope and hierarchy and so I have in some cases split up some of David's files to properly scope various function used in the vocoder. For example, the the three functions in the lsp{h,cpp} source files perform Line Spectrum Pairs are used in the main vocoder class and in the quantize functionality, so they have been separated.

Top level classes include:

- CCodec2: This class provides the primary interface for encoding and decoding. This class has a minimum count of public member functions that make up the external routines needed to be called to perform vocoding. There are a large number of private member functions that do the work. Nearly all other classes are are private members to this class and are used mostly to compartmentalize the various tasks.

- CKissFFT houses the forward and reverse Fourier transform code. Other than CCodec2, this is the only globally declared class.

- CQuantize, CNewamp1 and CNewamp2 are top level classes of a hiearchy of classes (CNewampbase and CQbase) that encapsulate all other procedural functions. The Newamp classes encapsulate the quantization functionality of the vocoder.

- Clpc and Cnlp are the Linear Prediction and Non-Linear Pitch functions and are declared as private members to the CCodec2 class.

There are several other significant changes to the Codec2 sources:

- Complex numbers are now implemented from the standard library (`#include <complex>`). One result is that complex math operators don't have to be specified, they are included in the standard library.

- All malloc-based memory allocation has been removed from the code. Array allocations are replaced with the vector standard container (`#include <vector>`) and code structures (`struct`) allocations are now delared as class member data, effectivey moving these objects from the heap to the stack.

- All tables of floating point numbers are contained in a single source file, codebooks.cpp.

## Building

Clone this repo and then see what it does:

```bash
git clone https://github.com/n7tae/c2-test.git
cd c2-test
make test
```

This will build the encoder and decoder and then run it with each of the supported Codec2 modes. The only operating system utility needed is `aplay`. Before and after the vocoding action, the original raw voice file will be played for comparison. When done, the coded files will be listed along with the raw voice file so you can see the size comparison. You can supply your own audio file if you want to see how it sounds. You can use `sox` to convert it to an 8000 Hz, 16-bit, monoral raw audio file needed by the encoder. See `sox --help`.
