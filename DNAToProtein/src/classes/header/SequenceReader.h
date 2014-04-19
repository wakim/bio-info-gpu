/*
 * MatrixReader.h
 *
 *  Created on: Jan 20, 2013
 *      Author: wakim
 */

#ifndef SEQUENCEREADER_H_
#define SEQUENCEREADER_H_

#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <exception>

namespace BioInfo {

	template <typename T>
	class SequenceReader {
	private:
		int fd;
		FILE *f;
		const char *path;
		bool closed;
	public:
		SequenceReader() : f(NULL), fd(-1), closed(true), path(NULL) {};
		SequenceReader(const char *_path) : f(NULL), fd(-1), closed(true), path(_path) {};
		virtual ~SequenceReader();
		bool open();
		T* read(int &size);
		T* raw_read(int& size);
		bool close();
		long get_file_size();
	};
};

#endif /* SEQUENCEREADER_H_ */
