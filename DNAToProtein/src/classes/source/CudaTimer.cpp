/*
 * CudaTimer.cpp
 *
 *  Created on: Jan 24, 2013
 *      Author: wakim
 */

#include "../header/CudaTimer.h"
#include "../../utils/cuda_utils.h"

CudaUtils::CudaTimer::CudaTimer() {
	CUDA_CHECK_RETURN(::cudaEventCreate(&_start));
	CUDA_CHECK_RETURN(::cudaEventCreate(&_stop));
	elapsed_time = -1;
}

CudaUtils::CudaTimer::~CudaTimer() {
	CUDA_CHECK_RETURN(cudaEventDestroy(_start));
	CUDA_CHECK_RETURN(cudaEventDestroy(_stop));
}

void CudaUtils::CudaTimer::start() {
	elapsed_time = -1;
	CUDA_CHECK_RETURN(cudaEventRecord(_start, DEFAULT_STREAM));
}

float CudaUtils::CudaTimer::stop() {
	if(elapsed_time != -1) {
		return elapsed_time;
	}

	CUDA_CHECK_RETURN(cudaEventRecord(_stop, DEFAULT_STREAM));
	CUDA_CHECK_RETURN(cudaEventSynchronize(_stop));

	CUDA_CHECK_RETURN(cudaEventElapsedTime(&elapsed_time, _start, _stop));

	return elapsed_time;
}
