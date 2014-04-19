/*
 * DNATranslatorImpl.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: wakim
 */

#include "../header/DNATranslatorImpl.h"
#include "../../utils/cuda_utils.h"

namespace BioInfo {

	template
	DNATranslatorImpl<char, 512>::~DNATranslatorImpl();

	template
	DNATranslatorImpl<char, 128>::~DNATranslatorImpl();

	template <typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::cleanup() {
		if(dna_seq_host != NULL) {
			delete dna_seq_host;
			dna_seq_host = NULL;
		}

		if(rna_seq_host != NULL) {
			delete rna_seq_host;
			rna_seq_host = NULL;
		}

		if(aa_seq_host != NULL) {
			delete aa_seq_host;
			aa_seq_host = NULL;
		}

		if(dna_seq_device != NULL) {
			CUDA_CHECK_RETURN(cudaFree(dna_seq_device));
			dna_seq_device = NULL;
		}

		if(rna_seq_device != NULL) {
			CUDA_CHECK_RETURN(cudaFree(rna_seq_device));
			rna_seq_device = NULL;
		}

		if(aa_seq_device != NULL) {
			CUDA_CHECK_RETURN(cudaFree(aa_seq_device));
			aa_seq_device = NULL;
		}
	}

	template <typename T, int max_threads>
	DNATranslatorImpl<T, max_threads>::~DNATranslatorImpl() {
		cleanup();
	}

	template <typename T, int max_threads>
	const int DNATranslatorImpl<T, max_threads>::get_max_chunk_size() const {
		cudaDeviceProp *prop = CudaUtils::get_cuda_device_properties();
		int result = max_threads * prop->maxGridSize[0] * 3;

		delete prop;

		return result;
	}

	template<typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::load_dna_seq_on_gpu() {
		CudaUtils::malloc_on_gpu<T>((void**) &dna_seq_device, multiple_3_size, debug);
		CUDA_CHECK_RETURN(cudaMemcpy(dna_seq_device, dna_seq_host, multiple_3_size * sizeof(T), cudaMemcpyHostToDevice));
	}

	template<typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::load_dna_seq_from_gpu() {
		if(dna_seq_device != NULL) {
			CUDA_CHECK_RETURN(cudaMemcpy(dna_seq_host, dna_seq_device, multiple_3_size * sizeof(T), cudaMemcpyDeviceToHost));
		}
	}

	template<typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::load_rna_seq_from_gpu() {
		if(rna_seq_device != NULL) {
			rna_seq_host = new T[multiple_3_size];
			CUDA_CHECK_RETURN(cudaMemcpy(rna_seq_host, rna_seq_device, size * sizeof(T), cudaMemcpyDeviceToHost));
		}
	}

	template<typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::load_aa_seq_from_gpu() {
		if(aa_seq_device != NULL) {
			aa_seq_host = new T[aa_seq_size];
			CUDA_CHECK_RETURN(cudaMemcpy(aa_seq_host, aa_seq_device, aa_seq_size * sizeof(T), cudaMemcpyDeviceToHost));
		}
	}

	template <typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::generate_rna_seq() {
		rna_seq_device = ::translation_wrapper_specialized<T, max_threads, false>(dna_seq_device, multiple_3_size, TranslationOperation::Operation::DNA_TO_RNA);
		//rna_seq_device = ::translation_wrapper<T, max_threads, false>(dna_seq_device, multiple_3_size, TranslationOperation::Operation::DNA_TO_RNA);
	}

	template <typename T, int max_threads>
	void DNATranslatorImpl<T, max_threads>::generate_aa_seq(bool complete) {

		if(! complete) {
			aa_seq_size = multiple_3_size / 3;
		} else {
			aa_seq_size = multiple_3_size;
		}

		if(rna_seq_device != NULL) {
			//aa_seq_device = ::translation_wrapper_specialized<T, max_threads, false>(rna_seq_device, multiple_3_size, TranslationOperation::Operation::RNA_TO_AA, complete);
			aa_seq_device = ::translation_wrapper<T, max_threads, false>(rna_seq_device, multiple_3_size, TranslationOperation::Operation::RNA_TO_AA, complete);
		} else {
			aa_seq_device = ::translation_wrapper_specialized<T, max_threads, false>(dna_seq_device, multiple_3_size, TranslationOperation::Operation::DNA_TO_AA, complete);
			//aa_seq_device = ::translation_wrapper<T, max_threads, false>(dna_seq_device, multiple_3_size, TranslationOperation::Operation::DNA_TO_AA, complete);
		}
	}

} /* namespace BioInfo */
