#ifndef FILEUTILS_H_
#define FILEUTILS_H_

#ifdef __linux__

#include <unistd.h>
#define error_to_string strerror
#define file_open open
#define file_read read
#define file_write write
#define file_close close

#elif _WIN32 | _WIN64

#include "unistd_win.h"
#include <share.h>
#define error_to_string strerror_s
#define file_open _sopen_s
#define file_read _read
#define file_write _write
#define file_close _close
#define WINDOWS 1

#endif

inline errno_t open_file(int &fd, const char* name, int flags, int shflag = _SH_DENYRW, int pmode = 0) {
	errno_t my_errno = 0;
	#ifdef WINDOWS
		errno = _sopen_s(&fd, name, flags, shflag, pmode);
		return errno;
	#else
		fd = open(name, flags);
		return errno;
	#endif
}
inline void print_error_string(errno_t the_errno) {
	char buffer[255];
	strerror_s(buffer, 255, errno);
	printf("%s\n", buffer);
}

#endif /* FILEUTILS_H_  */