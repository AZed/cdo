#ifndef  _CDIFORTRAN_H
#define  _CDIFORTRAN_H

/*******************************************************************************
 * Character buffer:
 */

#define CBUF_cfINT(N,A,B,X,Y,Z)		STRING_cfINT(N,A,B,X,Y,Z)
#define CBUF_cfSEP(T,  B)		STRING_cfSEP(T,B)
#define CBUF_cfN(  T,A)			STRING_cfN(T,A)
#define CBUF_cfSTR(N,T,A,B,C,D,E)	STRING_cfSTR(N,T,A,B,C,D,E)
#if defined(vmsFortran)
#   define CBUF_cfT(M,I,A,B,D)		A->dsc$a_pointer
#elif defined(CRAYFortran)
#   define CBUF_cfT(M,I,A,B,D)		_fcdtocp(A)
#else
#   define CBUF_cfT(M,I,A,B,D)		A
#endif

#endif
