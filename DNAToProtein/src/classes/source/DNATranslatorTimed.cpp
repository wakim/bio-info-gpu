/*
 * DNATranslatorTimed.cpp
 *
 *  Created on: Jan 24, 2013
 *      Author: wakim
 */

#include "../header/DNATranslatorTimed.h"

namespace BioInfo {

	template
	DNATranslatorTimed<char, 512>::~DNATranslatorTimed();

	template
	DNATranslatorTimed<char, 128>::~DNATranslatorTimed();

	template <typename T, int max_threads>
	DNATranslatorTimed<T, max_threads>::~DNATranslatorTimed() {
		delete instance;
	}

	template<typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::load_dna_seq_on_gpu() {
		timer.start();
		instance->load_dna_seq_on_gpu();
		float elapsed_time = timer.stop();
		times[TimerFunctions::LOAD_DNA_SEQ_ON_GPU] = elapsed_time;
	}

	template<typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::load_dna_seq_from_gpu() {
		timer.start();
		instance->load_dna_seq_from_gpu();
		float elapsed_time = timer.stop();
		times[TimerFunctions::LOAD_DNA_SEQ_FROM_GPU] = elapsed_time;
	}

	template<typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::load_rna_seq_from_gpu() {
		timer.start();
		instance->load_rna_seq_from_gpu();
		float elapsed_time = timer.stop();
		times[TimerFunctions::LOAD_RNA_SEQ_FROM_GPU] = elapsed_time;
	}

	template<typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::load_aa_seq_from_gpu() {
		timer.start();
		instance->load_aa_seq_from_gpu();
		float elapsed_time = timer.stop();
		times[TimerFunctions::LOAD_AA_SEQ_FROM_GPU] = elapsed_time;
	}

	template <typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::generate_rna_seq() {
		timer.start();
		instance->generate_rna_seq();
		float elapsed_time = timer.stop();
		times[TimerFunctions::GENERATE_RNA_SEQ] = elapsed_time;
	}

	template <typename T, int max_threads>
	void DNATranslatorTimed<T, max_threads>::generate_aa_seq(bool complete) {
		timer.start();
		instance->generate_aa_seq(complete);
		float elapsed_time = timer.stop();
		times[TimerFunctions::GENERATE_AA_SEQ] = elapsed_time;
	}

	template <typename T, int max_threads>
	float DNATranslatorTimed<T, max_threads>::get_sum_times() {
		float sum = 0.0f;

		printf("\n\n\nMetricas obtidas\n");
		for(int i = 0; i < TimerFunctions::LAST; ++i) {
			printf("%s - %f ms\n", TimerFunctions::TimerFunctionNames[i], times[i]);
			sum += times[i];
		}

		return sum;
	}
} /* namespace BioInfo */
