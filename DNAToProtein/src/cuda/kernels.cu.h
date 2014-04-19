/*
 * kernels.cu.h
 *
 *  Created on: Jan 21, 2013
 *      Author: wakim
 */

#ifndef KERNELS_CU_H_
#define KERNELS_CU_H_

#include "../utils/cuda_utils.h"

__device__ __constant__ const static char aa_code_device_lookup_table[4][4][4] =
{
	'K', 'N', 'N', 'K', // 0
	'I', 'I', 'I', 'M', // 1
	'T', 'T', 'T', 'T', // 2
	'R', 'S', 'S', 'R', // 3
	'*', 'Y', 'Y', '*', // 4
	'L', 'F', 'F', 'L', // 5
	'S', 'S', 'S', 'S', // 6
	'*', 'C', 'C', 'W', // 7
	'Q', 'H', 'H', 'Q', // 8
	'L', 'L', 'L', 'L', // 9
	'P', 'P', 'P', 'P', // 10
	'R', 'R', 'R', 'R', // 11
	'E', 'D', 'D', 'E', // 12
	'V', 'V', 'V', 'V', // 13
	'A', 'A', 'A', 'A', // 14
	'G', 'G', 'G', 'G', // 15
};


__device__ __constant__ const static char* aa_name_device_lookup_table[4][4][4] =
{
	"LYS", "ASN", "ASN", "LYS", // 0
	"ILE", "ILE", "ILE", "MET", // 1
	"THR", "THR", "THR", "THR", // 2
	"ARG", "SER", "SER", "ARG", // 3
	"***", "TYR", "TYR", "***", // 4
	"LEU", "PHE", "PHE", "LEU", // 5
	"SER", "SER", "SER", "SER", // 6
	"***", "CYS", "CYS", "TRP", // 7
	"GLN", "HIS", "HIS", "GLN", // 8
	"LEU", "LEU", "LEU", "LEU", // 9
	"PRO", "PRO", "PRO", "PRO", // 10
	"ARG", "ARG", "ARG", "ARG", // 11
	"GLU", "ASP", "ASP", "GLU", // 12
	"VAL", "VAL", "VAL", "VAL", // 13
	"ALA", "ALA", "ALA", "ALA", // 14
	"GLY", "GLY", "GLY", "GLY", // 15
};

namespace TranslationOperation {
	enum Operation {
		DNA_TO_RNA, DNA_TO_AA, RNA_TO_AA
	};
};

template <typename T, int max_threads, bool modify_input_sequence>
T* translation_wrapper(T* input_seq, int size, TranslationOperation::Operation op, bool complete = true);

template <typename T, int max_threads, bool modify_input_sequence>
T* translation_wrapper_specialized(T* input_seq, int size, TranslationOperation::Operation op, bool complete = true);

template <typename T>
__device__
int nucleotide_index(T nucleotide) {

	switch(nucleotide) {
		case 'A' : return 0;
		case 'U' :
		case 'T' : return 1;
		case 'C' : return 2;
		case 'G' : return 3;
	}

	return -1;
}

template <typename T>
__device__ void get_indexed_pos(T n1, T n2, T n3, int &p1, int &p2, int &p3) {
	p1 = nucleotide_index(n1);
	p2 = nucleotide_index(n2);
	p3 = nucleotide_index(n3);
}

template <typename T>
__device__
void apply_name_lookup_table(T* src, T* dest, int pos) {
	T nucleotide1 = src[pos];
	T nucleotide2 = src[pos + 1];
	T nucleotide3 = src[pos + 2];

	int lookup_index_1 = -1;
	int lookup_index_2 = -1;
	int lookup_index_3 = -1;

	get_indexed_pos(nucleotide1, nucleotide2, nucleotide3, lookup_index_1, lookup_index_2, lookup_index_3);

	// Pegar na lookup table o aa indexado pela cadeia de DNA
	if((lookup_index_1 & lookup_index_2 & lookup_index_3) != -1) {
		const char* aa_name = aa_name_device_lookup_table[lookup_index_1][lookup_index_2][lookup_index_3];

		dest[pos] = aa_name[0];
		dest[pos + 1] = aa_name[1];
		dest[pos + 2] = aa_name[2];
	} else {
		dest[pos] = '-';
		dest[pos + 1] = '-';
		dest[pos + 2] = '-';
	}
}

template <typename T>
__device__
void apply_code_lookup_table(T* src, T* dest, int pos) {
	int real_pos = pos / 3;
	T nucleotide1 = src[pos];
	T nucleotide2 = src[pos + 1];
	T nucleotide3 = src[pos + 2];

	int lookup_index_1 = -1;
	int lookup_index_2 = -1;
	int lookup_index_3 = -1;

	get_indexed_pos(nucleotide1, nucleotide2, nucleotide3, lookup_index_1, lookup_index_2, lookup_index_3);

	// Pegar na lookup table o aa indexado pela cadeia de DNA
	if((lookup_index_1 & lookup_index_2 & lookup_index_3) != -1) {
		dest[real_pos] = aa_code_device_lookup_table[lookup_index_1][lookup_index_2][lookup_index_3];
	} else {
		dest[real_pos] = '-';
	}
}

class DNARToRNAPositionCalculator {
public:
	inline __device__ int operator() (unsigned int tx, unsigned int bx, unsigned int bdx, int size, bool &in_bounds) {
		int position = (bx * bdx + tx) * 3;
		
		in_bounds = position < size;
		return position;
	}
};

class XNAToProteinPositionCalculator {
public:
	inline __device__ int operator() (unsigned int tx, unsigned int bx, unsigned int bdx, int size, bool &in_bounds) {
		int position = (bx * bdx + tx) * 3;

		in_bounds = position < size;

		return position;
	}
};

template <typename T>
class DNAToProteinTranslationOperator {
public:
	DNAToProteinTranslationOperator(bool complete) { this->complete = complete; }
	inline __device__ __host__ void operator() (T* src, int src_size, T* dest, int dest_size, int pos) {
		if(complete) {
			apply_name_lookup_table(src, dest, pos);
		} else {
			apply_code_lookup_table(src, dest, pos);
		}
	}
private:
	 bool complete;
};

template <typename T>
class DNAToRNATranslationOperator {
public:
	inline __device__ void operator() (T* src, int src_size, T* dest, int dest_size, int pos) {
		T nucleotide1 = src[pos] & 223;
		T nucleotide2 = src[pos+1] & 223;
		T nucleotide3 = src[pos+2] & 223;

		if(nucleotide1 == 'T') {
			nucleotide1 = 'U';
		}

		if(nucleotide2 == 'T') {
			nucleotide2 = 'U';
		}

		if(nucleotide3 == 'T') {
			nucleotide3 = 'U';
		}

		// upper-case
		dest[pos] = nucleotide1;
		dest[pos+1] = nucleotide2;
		dest[pos+2] = nucleotide3 ;
	}
};

template <typename T>
class RNAToProteinTranslationOperator {
public:
	RNAToProteinTranslationOperator(bool complete) { this->complete = complete; }
	inline __device__ __host__ void operator() (T* src, int src_size, T* dest, int dest_size, int pos) {
		if(complete) {
			apply_name_lookup_table(src, dest, pos);
		} else {
			apply_code_lookup_table(src, dest, pos);
		}
	}
private:
	bool complete;
};

#endif /* KERNELS_CU_H_ */
