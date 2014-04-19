/*
 * cuda_utils.h
 *
 *  Created on: Jan 22, 2013
 *      Author: wakim
 */

#ifndef CUDA_UTILS_H_
#define CUDA_UTILS_H_

#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <typeinfo>

#define CUDA_CHECK_RETURN(value) {								\
	cudaError_t _m_cudaStat = value;							\
	if (_m_cudaStat != cudaSuccess) {							\
		fprintf(stderr, "Error %s at line %d in file %s\n",		\
		cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
		exit(1);												\
	}															\
}

namespace CudaUtils {
	inline unsigned int next_pow2(unsigned int x) {
		--x;
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;
		return ++x;
	}

	inline unsigned int best_sqrt(int x) {
		float f_x = static_cast<float>(x);

		f_x = sqrt(f_x);
		x = static_cast<int>(ceil(f_x));

		return x;
	}

	inline void compute_num_threads_blocks(dim3& grid, dim3& block, int max_threads, int size, bool check = true) {

		int sqrt_max_threads = max_threads;

		if(sqrt_max_threads > size) {
			block.x = size;
		} else {
			block.x = sqrt_max_threads;
		}

		grid.x = static_cast<unsigned int>(size / block.x);

		grid.x += (size % block.x == 0 ? 0 : 1);

		if(check) {
			cudaDeviceProp prop;
			int device;

			CUDA_CHECK_RETURN(cudaGetDevice(&device));
			CUDA_CHECK_RETURN(cudaGetDeviceProperties(&prop, device));

			int usage = block.x * grid.x;
			int capacity = prop.maxThreadsDim[0] * prop.maxGridSize[0];
			
			if(usage > capacity) {
				printf("Capacidade de processamento excedida, entre com um sequencia menor !\n");
				exit(-1);
			}

			double usage_d = static_cast<double>(usage) / static_cast<double>(capacity);
			printf("Uso de processamento %lf%%\n", usage_d);
		}
	}

	inline cudaDeviceProp* get_cuda_device_properties() {
		cudaDeviceProp *prop = new cudaDeviceProp();
		int device;

		CUDA_CHECK_RETURN(cudaGetDevice(&device));
		CUDA_CHECK_RETURN(cudaGetDeviceProperties(prop, device));

		return prop;
	}

	template <typename T>
	T* generate_reverse_linear_matrix(unsigned int rows, unsigned int cols) {
		T* matrix = new T[rows * cols];
		int max = rows * cols + 1;

		for(unsigned int i = 0; i < rows; ++i) {
			for(unsigned int j = 0; j < cols; ++j) {
				matrix[i * cols + j] = static_cast<T>(max--);
			}
		}

		return matrix;
	}

	template <typename T>
	const char* get_format() {
		if(typeid(T) == typeid(int)) {
			return "%d";
		} else if(typeid(T) == typeid(float)) {
			return "%1.2f";
		} else if(typeid(T) == typeid(double)) {
			return "%1.2lf";
		} else if(typeid(T) == typeid(char)) {
			return "%c";
		}
		return "";
	}

	template <typename T, bool line_number>
	void load_and_print(void* matrix_device, int rows, int cols) {
		int size = rows * cols;
		int size_bytes = size * sizeof(T);
		T* matrix_host = new T[size];
		const char* format = get_format<T>();

		if(line_number) {
			printf("format %s\n", format);
		}

		CUDA_CHECK_RETURN(cudaMemcpy(matrix_host, matrix_device, size_bytes, cudaMemcpyDeviceToHost));

		for(int i = 0; i < rows; ++i) {
			if(line_number) {
				printf("%d)", i);
			}
			for(int j = 0; j < cols; ++j) {
				printf(format, matrix_host[i * cols + j]);
			}
			printf("\n");
		}

		delete matrix_host;
	}

	template <typename T>
	void malloc_on_gpu(void** device_ptr, int size, bool debug = false) {
		size_t free, total;

		if(debug) {
			CUDA_CHECK_RETURN(cudaMemGetInfo(&free, &total));
			printf("Informacao sobre memoria (Total: %u, Livre: %u)\n", total, free);
		}

		CUDA_CHECK_RETURN(cudaMalloc(device_ptr, size * sizeof(T)));

		if(debug) {
			CUDA_CHECK_RETURN(cudaMemGetInfo(&free, &total));
			printf("Informacao sobre memoria (Total: %u, Livre: %u)\n", total, free);
		}

		if(debug) {
			double memory_usage = static_cast<double>(free) / static_cast<double>(total);
			printf("Usando %1.lf%% do total de memoria disponivel na GPU\n", memory_usage);
		}
	}
}

#endif /* CUDA_UTILS_H_ */
