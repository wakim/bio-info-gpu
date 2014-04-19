/*
 * MatrixReader.h
 *
 *  Created on: Jan 20, 2013
 *      Author: wakim
 */

#ifndef SEQUENCEWRITER_H_
#define SEQUENCEWRITER_H_

#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <exception>

namespace BioInfo {
	template <class T>
	class SequenceWriter {
	private:
		int fd;
		FILE *f;
		bool closed;
		const char *path;
	public:
		SequenceWriter() : f(NULL), fd(-1), path(NULL), closed(true){};
		SequenceWriter(const char *_path) : f(NULL), fd(-1), path(_path), closed(true) {};
		virtual ~SequenceWriter();
		bool open();
		bool write(const T* stream, const int size);
		bool close();
	};
};

#endif /* SEQUENCEWRITER_H_ */
