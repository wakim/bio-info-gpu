#include "../header/SequenceWriter.h"

namespace BioInfo {

	template
	SequenceWriter<char>::~SequenceWriter();

	template
	bool SequenceWriter<char>::open();

	template
	bool SequenceWriter<char>::write(const char* stream, const int size);

	template
	bool SequenceWriter<char>::close();

	template <typename T>
	SequenceWriter<T>::~SequenceWriter() {
		close();
	}

	template <typename T>
	bool SequenceWriter<T>::open() {
		fopen_s(&f, path, "w");
		
		if(f == NULL || errno != 0) {
			fprintf(stderr, "Erro ao abrir arquivo %s errno %d\n", path, errno);
			closed = true;
			return false;
		}
		
		fd = _fileno(f);
		closed = false;

		return true;
	}

	template <typename T>
	bool SequenceWriter<T>::close() {
		if(closed) {
			return true;
		}

		if(fclose(f) != 0) {
			fprintf(stderr, "Erro ao fechar arquivo %s\n", path);
			return false;
		}

		closed = true;

		return true;
	}

	template<typename T>
	bool SequenceWriter<T>::write(const T* stream, const int size) {
		int size_bytes = size * sizeof(T);

		if(closed) {
			return false;
		}

		if(fwrite(stream, sizeof(T), size, f) != size) {
			fprintf(stderr, "Erro ao escrever a sequencia no arquivo %s\n", path);
			return false;
		}

		return true;
	}
}
