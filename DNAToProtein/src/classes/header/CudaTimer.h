/*
 * CudaTimer.h
 *
 *  Created on: Jan 24, 2013
 *      Author: wakim
 */

#ifndef CUDATIMER_H_
#define CUDATIMER_H_

#include <cuda_runtime.h>

namespace CudaUtils {

	static const cudaStream_t DEFAULT_STREAM = 0x00;

	class CudaTimer {
		public:
			CudaTimer();
			virtual ~CudaTimer();
			void start();
			float stop();
		private:
			cudaEvent_t _start;
			cudaEvent_t _stop;
			float elapsed_time;
	};
}

#endif /* CUDATIMER_H_ */
