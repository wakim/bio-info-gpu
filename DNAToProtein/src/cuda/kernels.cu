/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include <stdio.h>
#include <typeinfo>
#include "shared_memory.cu.h"
#include "kernels.cu.h"

template <typename T>
__global__ void dna_to_rna_sequence_translator_kernel(T* input_seq, int input_size, T* output_seq, int output_size) {

	unsigned int tx = threadIdx.x;
	unsigned int bx = blockIdx.x;
	unsigned int bdx = blockDim.x;
	int position = (bx * bdx + tx) * 3;
	bool in_bounds = position < input_size;

	if(in_bounds) {
		T nucleotide1 = input_seq[position];
		T nucleotide2 = input_seq[position+1];
		T nucleotide3 = input_seq[position+2];

		if(nucleotide1 == 'T' || nucleotide1 == 't') {
			nucleotide1 = 'U';
		}

		if(nucleotide2 == 'T' || nucleotide2 == 't') {
			nucleotide2 = 'U';
		}

		if(nucleotide3 == 'T' || nucleotide3 == 't') {
			nucleotide3 = 'U';
		}

		// upper-case
		output_seq[position] = nucleotide1 &= 223;
		output_seq[position+1] = nucleotide2 &= 223;
		output_seq[position+2] = nucleotide3 &= 223;
	}
}

template <typename T>
__global__ void dna_to_aa_sequence_translator_kernel(T* input_seq, int input_size, T* output_seq, int output_size, bool complete) {

	unsigned int tx = threadIdx.x;
	unsigned int bx = blockIdx.x;
	unsigned int bdx = blockDim.x;
	int position = (bx * bdx + tx) * 3;
	bool in_bounds = position < input_size;

	if(in_bounds) {
		if(complete) {
			apply_name_lookup_table(input_seq, output_seq, position);
		} else {
			apply_code_lookup_table(input_seq, output_seq, position);
		}
	}
}

template <typename T>
__global__ void rna_to_aa_sequence_translator_kernel(T* input_seq, int input_size, T* output_seq, int output_size, bool complete) {

	unsigned int tx = threadIdx.x;
	unsigned int bx = blockIdx.x;
	unsigned int bdx = blockDim.x;
	int position = (bx * bdx + tx) * 3;
	bool in_bounds = position < input_size;

	if(in_bounds) {
		if(complete) {
			apply_name_lookup_table(input_seq, output_seq, position);
		} else {
			apply_code_lookup_table(input_seq, output_seq, position);
		}
	}
}

template <typename T, typename TranslationOperator, typename PositionCalculator>
__global__ void general_sequence_translator_kernel(T* input_seq, int input_size, T* output_seq, int output_size, TranslationOperator op1, PositionCalculator op2) {

	unsigned int tx = threadIdx.x;
	unsigned int bx = blockIdx.x;
	unsigned int bdx = blockDim.x;
	bool in_bounds = false;
	unsigned int pos = op2(tx, bx, bdx, input_size, in_bounds);

	if(in_bounds) {
		op1(input_seq, input_size, output_seq, output_size, pos);
	}
}

template <typename T, int max_threads, bool modify_input_sequence, typename TranslationOperator, typename PositionCalculator>
T* _translation_wrapper(T* input_seq, int input_size, int output_size, TranslationOperator op1, PositionCalculator op2) {

	dim3 grid, block;
	T *output_seq = input_seq;
	
	CudaUtils::compute_num_threads_blocks(grid, block, max_threads, output_size, false);

	if(! modify_input_sequence) {
		output_seq = NULL;
		CudaUtils::malloc_on_gpu<T>((void**) &output_seq, input_size);
	}

	general_sequence_translator_kernel<T, TranslationOperator, PositionCalculator>
		<<<grid, block>>>(input_seq, input_size, output_seq, output_size, op1, op2);

	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	//printf("\n\noutput_seq\n\n");
	//CudaUtils::load_and_print<char, true>(output_seq, 1, input_size);

	return output_seq;
}

template <typename T, int max_threads, bool modify_input_sequence>
T* translation_wrapper(T* input_seq, int size, TranslationOperation::Operation op, bool complete) {

	int size_div_3 = size / 3;

	switch(op) {
		case TranslationOperation::DNA_TO_RNA:
			return _translation_wrapper<T, max_threads, modify_input_sequence, DNAToRNATranslationOperator<T>, DNARToRNAPositionCalculator>
				(input_seq, size, size_div_3, DNAToRNATranslationOperator<T>(), DNARToRNAPositionCalculator());
		case TranslationOperation::DNA_TO_AA:
			return _translation_wrapper<T, max_threads, modify_input_sequence, DNAToProteinTranslationOperator<T>, XNAToProteinPositionCalculator>
				(input_seq, size, size_div_3, DNAToProteinTranslationOperator<T>(complete), XNAToProteinPositionCalculator());
		case TranslationOperation::RNA_TO_AA:
			return _translation_wrapper<T, max_threads, modify_input_sequence, RNAToProteinTranslationOperator<T>, XNAToProteinPositionCalculator>
				(input_seq, size, size_div_3, RNAToProteinTranslationOperator<T>(complete), XNAToProteinPositionCalculator());
	}

	return NULL;
}

template <typename T, int max_threads, bool modify_input_sequence>
T* translation_wrapper_specialized(T* input_seq, int size, TranslationOperation::Operation op, bool complete) {
	int output_size = size / 3;
	dim3 grid, block;
	T *output_seq = input_seq;
	
	CudaUtils::compute_num_threads_blocks(grid, block, max_threads, output_size, true);

	if(! modify_input_sequence) {
		output_seq = NULL;
		CudaUtils::malloc_on_gpu<T>((void**) &output_seq, output_size);
	}

	switch(op) {
		case TranslationOperation::DNA_TO_RNA:
			dna_to_rna_sequence_translator_kernel<T> <<<grid, block>>>(input_seq, size, output_seq, output_size);
			break;
		case TranslationOperation::DNA_TO_AA:
			dna_to_aa_sequence_translator_kernel<T> <<<grid, block>>>(input_seq, size, output_seq, output_size, complete);
			break;
		case TranslationOperation::RNA_TO_AA:
			rna_to_aa_sequence_translator_kernel<T> <<<grid, block>>>(input_seq, size, output_seq, output_size, complete);
			break;
	}
	
	cudaError_t last_error = cudaDeviceSynchronize();
	CUDA_CHECK_RETURN(last_error);

	//printf("\n\noutput_seq\n\n");
	//CudaUtils::load_and_print<char, true>(output_seq, 1, input_size);

	return output_seq;
}

// 512 threads

template
char* translation_wrapper_specialized<char, 512, true> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper_specialized<char, 512, false> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper<char, 512, true> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper<char, 512, false> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

// 128 threads

template
char* translation_wrapper_specialized<char, 128, true> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper_specialized<char, 128, false> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper<char, 128, true> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);

template
char* translation_wrapper<char, 128, false> (char* input_seq, int size, TranslationOperation::Operation op, bool complete);