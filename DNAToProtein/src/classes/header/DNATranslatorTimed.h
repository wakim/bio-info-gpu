/*
 * DNATranslatorTimed.h
 *
 *  Created on: Jan 24, 2013
 *      Author: wakim
 */

#ifndef DNATRANSLATORTIMED_H_
#define DNATRANSLATORTIMED_H_

#include "DNATranslator.h"
#include "DNATranslatorImpl.h"
#include "CudaTimer.h"

namespace BioInfo {

	namespace TimerFunctions {
		enum {
			LOAD_DNA_SEQ_ON_GPU,
			LOAD_DNA_SEQ_FROM_GPU, LOAD_RNA_SEQ_FROM_GPU, LOAD_AA_SEQ_FROM_GPU,
			GENERATE_RNA_SEQ, GENERATE_AA_SEQ,

			LAST
		};

		static const char* TimerFunctionNames[] = {
			"LOAD_DNA_SEQ_ON_GPU",
			"LOAD_DNA_SEQ_FROM_GPU", "LOAD_RNA_SEQ_FROM_GPU", "LOAD_AA_SEQ_FROM_GPU",
			"GENERATE_RNA_SEQ", "GENERATE_AA_SEQ"
		};
	}

	/**
	 * Decorator para medir o tempo de execução de cada método.
	 * O certo é ter uma interface mais alto nível que DNATranslator
	 */
	template <typename T, int max_threads>
	class DNATranslatorTimed : public DNATranslator<T, max_threads> {
		public:
			DNATranslatorTimed(DNATranslatorImpl<T, max_threads>* _instance) : instance(_instance){
				memset(&times, 0, TimerFunctions::LAST * sizeof(float));
			}
			virtual ~DNATranslatorTimed();
			
			virtual inline void set_dna_seq(T* dna_seq, int size) {
				instance->set_dna_seq(dna_seq, size);
			}

			virtual inline const T* get_dna_seq() const { return instance->get_dna_seq(); }
			virtual inline const T* get_rna_seq() const { return instance->get_rna_seq(); }
			virtual inline const T* get_aa_seq(int &size) const { return instance->get_aa_seq(size); }
			virtual const int get_max_chunk_size() const { return instance->get_max_chunk_size(); }

			virtual void load_dna_seq_on_gpu();
			virtual void load_dna_seq_from_gpu();
			virtual void load_rna_seq_from_gpu();
			virtual void load_aa_seq_from_gpu();

			virtual void generate_rna_seq();
			virtual void generate_aa_seq(bool complete);

			virtual inline void set_debug(bool debug) { instance->set_debug(debug); }

			virtual float get_sum_times();
		private:
			DNATranslatorImpl<T, max_threads>* instance;
			CudaUtils::CudaTimer timer;
			float times[TimerFunctions::LAST];
	};

} /* namespace BioInfo */
#endif /* DNATRANSLATORTIMED_H_ */
