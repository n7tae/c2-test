# c2-cl-ser

## Description

This is a C++ implementation of David Row's Codec2. Many changes were performed to translate David's C code to be supported in a C++11 standard. It is interesting to note that a significant portion of the Codec2's source code are tables of floating point numbers. These numbers are used to process the input and output of forward and reverse Fourier Transforms[^1]. All operational code has been encapsulated in a few classes. The classes selected are based on functionality and so I have in some cases split up some of David's files to bring functionality. For example, the the three functions in the lsp{h,cpp} source files perform Line Spectrum Pairs are used in the top level class and in the quantize functionality, so they have been separated.

Top level classes include:

- CCodec2: This class provides the primary interface for encoding and decoding. This class has a minimum count of public members that make up the external routines needed to be called to perform vocoding. There are a number of private members that do the work. All other classes are are private members to this class and are used mostly to compartmentalize the various tasks.

- Clpc and Cnlp: The Linear Prediction and Non-Linear Pitch functions are classes declared as private members to the CCodec2 class.

- CQuantize, CNewamp1 and CNewamp2 are top level classes of a hiearchy of classes (CNewampbase and CQbase) that encapsulate all other procedural functions.

- All tables of floating point numbers are contained in a single source file, codebooks.cpp.

[^1]: [Fourier Transform](https://en.wikipedia.org/wiki/Fourier_transform) is used to translate time-series data (like a sound wave) into frequency-series data (like a waterfall display) and *vis-versa*.
