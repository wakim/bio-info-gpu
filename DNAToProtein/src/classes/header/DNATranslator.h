
#ifndef DNATRANSLATOR_H_
#define DNATRANSLATOR_H_

namespace BioInfo {
	template <typename T, int max_threads>
	class DNATranslator {
	public:
		virtual ~DNATranslator() {};
		virtual void set_dna_seq(T* dna_seq, int size) = 0;
		virtual const T* get_dna_seq() const = 0;
		virtual const T* get_rna_seq() const = 0;
		virtual const T* get_aa_seq(int &size) const = 0;
		virtual const int get_max_chunk_size() const = 0;

		/**
		 * Carrega a sequencia de DNA do Host para o Device
		 */
		virtual void load_dna_seq_on_gpu() = 0;

		/**
		 * Carrega a sequencia de DNA do Device para o Host
		 */
		virtual void load_dna_seq_from_gpu() = 0;

		/**
		 * Carrega a sequencia de RNA do Device para o Host
		 */
		virtual void load_rna_seq_from_gpu() = 0;

		/**
		 * Carrega a sequencia de Aminoácidos do Device para o Host
		 */
		virtual void load_aa_seq_from_gpu() = 0;

		/**
		 * Transforma uma sequência de DNA em RNA na GPU.
		 */
		virtual void generate_rna_seq() = 0;

		/**
		 * Transforma uma sequência de RNA em AA na GPU.
		 */
		virtual void generate_aa_seq(bool complete) = 0;

		virtual void set_debug(bool debug) = 0;
	};
}
#endif /* DNATRANSLATOR_H_ */
