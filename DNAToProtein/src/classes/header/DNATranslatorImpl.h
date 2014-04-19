/*
 * DNATranslatorImpl.h
 *
 *  Created on: Jan 21, 2013
 *      Author: wakim
 */

#ifndef DNATRANSLATORIMPL_H_
#define DNATRANSLATORIMPL_H_

#include <stdio.h>
#include <exception>
#include "DNATranslator.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "../../cuda/kernels.cu.h"

namespace BioInfo {

	template <typename T, int max_threads>
	class DNATranslatorImpl : public DNATranslator<T, max_threads> {
	private:
		T* dna_seq_host, *dna_seq_device;
		T* rna_seq_host, *rna_seq_device;
		T* aa_seq_host, *aa_seq_device;

		int size, multiple_3_size, aa_seq_size;
		dim3 block, grid;
		bool complete; // Sequencia completa ou apenas uma letra
		bool debug;
		void cleanup();
	public:
		DNATranslatorImpl(T* _dna_seq_host) : dna_seq_device(NULL), dna_seq_host(_dna_seq_host),
			rna_seq_host(NULL), rna_seq_device(NULL), aa_seq_host(NULL), aa_seq_device(NULL) {}
		DNATranslatorImpl() : dna_seq_device(NULL), dna_seq_host(NULL),
			rna_seq_host(NULL), rna_seq_device(NULL), aa_seq_host(NULL), aa_seq_device(NULL) {}

		virtual ~DNATranslatorImpl();

		virtual inline void set_dna_seq(T* dna_seq, int size) {
			int remainder = (size % 3);
			
			cleanup();

			dna_seq_host = dna_seq;

			this->multiple_3_size = size - remainder;
			this->size = size;
		}
		virtual inline const T* get_dna_seq() const {return dna_seq_host; }
		virtual inline const T* get_rna_seq() const { return rna_seq_host; }
		virtual inline const T* get_aa_seq(int &size) const { size = aa_seq_size; return aa_seq_host; }
		virtual const int get_max_chunk_size() const;

		virtual void load_dna_seq_on_gpu();
		virtual void load_dna_seq_from_gpu();
		virtual void load_rna_seq_from_gpu();
		virtual void load_aa_seq_from_gpu();

		virtual void generate_rna_seq();
		virtual void generate_aa_seq(bool complete);

		virtual inline void set_debug(bool debug) { this->debug = debug; }
	};

} /* namespace BioInfo */
#endif /* DNATRANSLATORIMPL_H_ */
