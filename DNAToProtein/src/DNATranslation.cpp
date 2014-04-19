#include <cstdio>
#include <cstring>
#include <ctime>

#include "classes/header/SequenceReader.h"
#include "classes/header/SequenceWriter.h"
#include "classes/header/DNATranslator.h"
#include "classes/header/DNATranslatorImpl.h"
#include "classes/header/DNATranslatorTimed.h"

#if defined _WIN32 || defined _WIN64
#define WAIT puts("Aperte ENTER para sair..."); getchar();
#else
#define WAIT puts("Aperte ENTER para sair...");
#endif

void print_complete_sequence(const char* seq, int size, const int start_index);
void print_compact_sequence(const char* seq, int size, const int start_index);
int get_number_of_chunks(const long sequence_size, int& max_chunk_size);

void print_usage(char* program_name) {
	printf("Uso: %s arquivo_entrada [tipo] [arquivo_saida] [debug] (tipo := completo | compacto)\n", program_name);
	printf("No tipo completo cada aa eh representado com tres letras (e.g: MET), no modo compacto apenas uma (e.g: M)\n");
}

int main(int argc, char **argv) {
	char* sequence, *input_path = NULL, *type = "completo";
	char output_path[255];
	bool debug_flag = false, complete = true;
	const int max_threads = 128;
	clock_t start, stop;

	start = clock();
	
	if(argc < 2) {
		print_usage(argv[0]);
		exit(1);
	}
	
	input_path = argv[1];

	if(argc > 2) {
		type = argv[2];
	}

	char max_threads_buffer[4];

	itoa(max_threads, max_threads_buffer, 10);

	if(argc > 3) {
		strcpy_s<255>(output_path, max_threads_buffer);
		strcat_s<255>(output_path, argv[3]);
	} else {
		strcpy_s<255>(output_path, input_path);
		strcat_s<255>(output_path, max_threads_buffer);
		strcat_s<255>(output_path, ".out");
	}

	if(argc > 4) {
		debug_flag = strlen(argv[4]) == 5 && strncmp("debug", argv[4], 5) == 0;
	}

	BioInfo::DNATranslator<char, max_threads>* translator = new BioInfo::DNATranslatorTimed<char, max_threads>(new BioInfo::DNATranslatorImpl<char, max_threads>());
	BioInfo::DNATranslatorTimed<char, max_threads> *real_translator = static_cast<BioInfo::DNATranslatorTimed<char, max_threads>*>(translator);
	BioInfo::SequenceReader<char> reader(input_path);
	BioInfo::SequenceWriter<char> writer(output_path);

	if(strlen(type) == 8) {
		if(strncmp("completo", type, 8) == 0) {
			complete = true;
		} else if(strncmp("compacto", type, 8) == 0) {
			complete = false;
		}
	}

	if(! reader.open()) {
		delete translator;
		exit(EXIT_FAILURE);
	}

	long sequence_size = reader.get_file_size();
	int remainder = sequence_size % 3;
	
	if(remainder) {
		printf("Alerta, a cadeia de entrada nao eh multipla de 3. Serao ignoradas as %d ultimas bases\n", remainder);
	}
	
	// Calcula quantas vezes devemos ir para a GPU processar a sequencia inteira.
	int max_chunk_size = translator->get_max_chunk_size();
	int chunk_size = max_chunk_size;
	int number_of_chunks = get_number_of_chunks(sequence_size, chunk_size);
	long offset = 0, total_processed = 0;

	if(number_of_chunks > 1) {
		printf("A sequencia sera processada %d vezes em pedacos de %d bytes\n", number_of_chunks, max_chunk_size);
	} else {
		printf("A sequencia tem %d bytes, o maior pedaco que pode ser tratado de uma vez e %d\n", chunk_size, max_chunk_size);
	}

	writer.open();

	for(int i = 0; i < number_of_chunks; ++i) {
		sequence = reader.raw_read(chunk_size);

		offset += chunk_size;
		total_processed += chunk_size;

		translator->set_debug(debug_flag);
		translator->set_dna_seq(sequence, chunk_size);
		translator->load_dna_seq_on_gpu();
		translator->generate_aa_seq(complete);
		translator->load_aa_seq_from_gpu();

		int aa_seq_size = chunk_size;

		const char* aa_seq = translator->get_aa_seq(aa_seq_size);

		if(! writer.write(aa_seq, aa_seq_size)) {
			fprintf(stderr, "Erro ao escrever saida em %s\n", output_path);
			break;
		}

		printf("Total Processado %ld/%ld pb\n", total_processed, sequence_size);
	}

	if(! reader.close()) {
		exit(EXIT_FAILURE);
	}

	if(! writer.close()) {
		exit(EXIT_FAILURE);
	}

	printf("Arquivo %s gerado\n", output_path);

	if(debug_flag) {
		float time = real_translator->get_sum_times();
		printf("Tempo total gasto na GPU: %f ms\n", time);
	}

	delete translator;
	
	stop = clock();

	double diff_time = ((double)(stop - start)) / CLOCKS_PER_SEC;

	printf("Tempo total gasto na cpu: %lf ms (%lf s)\n", diff_time * 1000.0, diff_time);

	WAIT
	return 0;
}

void print_complete_sequence(const char* seq, int size, const int start_index) {
	const char *p_seq = seq + start_index;
	printf("%s\n", p_seq);
}

void print_compact_sequence(const char* seq, int size, const int start_index) {
	const char *p_seq = seq + start_index;
	printf("%s\n", p_seq);
}

int get_number_of_chunks(const long sequence_size, int& max_chunk_size) {
	long l_max_chunk_size = static_cast<long>(max_chunk_size);
	
	if(sequence_size <= l_max_chunk_size) {
		max_chunk_size = static_cast<int>(sequence_size);
		return 1;
	}

	long remainder = l_max_chunk_size % 3;
	double f_chunk_size = static_cast<double>(max_chunk_size - remainder);
	double f_sequence_size = static_cast<double>(sequence_size);
	double f_number_of_chunks;

	f_number_of_chunks = ceil(f_sequence_size / f_chunk_size);

	return static_cast<int>(f_number_of_chunks);
}