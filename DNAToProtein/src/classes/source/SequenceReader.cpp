#include "../header/SequenceReader.h"

namespace BioInfo {

	template
	SequenceReader<char>::~SequenceReader();

	template
	bool SequenceReader<char>::open();

	template
	char* SequenceReader<char>::read(int &size);

	template
	char* SequenceReader<char>::raw_read(int& size);

	template
	bool SequenceReader<char>::close();

	template <typename T>
	SequenceReader<T>::~SequenceReader() {
		close();
	}

	template <typename T>
	bool SequenceReader<T>::open() {
		fopen_s(&f, path, "r");

		if(f == NULL || errno != 0) {
			fprintf(stderr, "Erro ao abrir arquivo %s errno %d\n", path, errno);
			return false;
		}

		fd = _fileno(f);
		closed = false;
		return true;
	}

	template <typename T>
	bool SequenceReader<T>::close() {
		if(closed) {
			return true;
		}

		if(fclose(f) != 0) {
			fprintf(stderr, "Erro ao fechar o arquivo %s\n", path);
			closed = false;
			return false;
		}

		closed = true;

		return true;
	}

	template<typename T>
	long SequenceReader<T>::get_file_size() {
		struct stat file_stats;

		if(closed) {
			return NULL;
		}

		if(fstat(fd, &file_stats) != 0) {
			fprintf(stderr, "Erro ao recuperar tamanho do arquivo %s\n", path);
			return -1;
		}

		return file_stats.st_size;
	}

	template<typename T>
	T* SequenceReader<T>::raw_read(int& size) {
		T* result = NULL;

		try {
			result = new T[size + 1];
		} catch(std::exception& e) {
			fprintf(stderr, "Erro ao alocar memoria para a sequencia %s\n", e.what());
			return NULL;
		}

		size = fread(result, sizeof(T), size, f);

		result[size] = 0;
		return result;
	}

	template<typename T>
	T* SequenceReader<T>::read(int& size) {
		size = get_file_size();

		if((size) <= 0) {
			fprintf(stderr, "Erro de leitura sequence_size <= 0\n");
			return NULL;
		}

		T* result = NULL;

		try {
			result = new T[size + 1];
		} catch(std::exception& e) {
			fprintf(stderr, "Erro ao alocar memoria para a sequencia %s\n", e.what());
			return NULL;
		}

		fread(result, sizeof(T), size, f);

		result[size] = 0;
		return result;
	}
}
