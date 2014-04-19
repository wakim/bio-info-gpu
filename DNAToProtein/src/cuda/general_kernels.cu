
#include "shared_memory.cu.h"

// Copyright http://www.evl.uic.edu/aej/525/
// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
template <typename T>
__global__ void transpose(T *odata, T *idata, int width, int height) {
	
	SharedMemory<T> shared_memory;
	T* block = shared_memory.getPointer();
	unsigned int BLOCK_DIM = blockDim.x;
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;

	if((xIndex < width) && (yIndex < height)) {
		unsigned int index_in = yIndex * width + xIndex;
		unsigned int index_block = threadIdx.y * BLOCK_DIM + threadIdx.x;

		block[index_block] = idata[index_in];
	}

	__syncthreads();

	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
	yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;

	if((xIndex < height) && (yIndex < width)) {
		unsigned int index_out = yIndex * height + xIndex;
		unsigned int index_block = threadIdx.x * BLOCK_DIM + threadIdx.y;

		odata[index_out] = block[index_block];
	}
}

template <typename T, typename N, typename NormalizeOperation>
__global__ void multiple_normalize(N* output_matrix, T* input_matrix, T* max_block_matrix, int rows, int cols,
	int block_matrix_cols, NormalizeOperation op) {

	unsigned int tx = threadIdx.x, ty = threadIdx.y;
	unsigned int bx = blockIdx.x, by = blockIdx.y;
	unsigned int bdx = blockDim.x, bdy = blockDim.y;
	unsigned int col = (bx * bdx) + tx;
	unsigned int row = (by * bdy) + ty;
	unsigned int matrix_pos = row * cols + col;
	bool in_bounds = row < rows && col < cols;

	if(in_bounds) {
		unsigned int block_matrix_pos = (row * block_matrix_cols) + (block_matrix_cols - 1);
		op(output_matrix, input_matrix, max_block_matrix[block_matrix_pos], matrix_pos);
	}
}

template <typename T, typename N, typename TFIDFOperation>
__global__ void transform_tf_idf(N* output_matrix, T* input_matrix, int rows, int cols,
	T* dfi_matrix, T* max_fj_matrix, int dfi_matrix_cols, int max_fj_matrix_cols,
	TFIDFOperation op) {

	unsigned int tx = threadIdx.x, ty = threadIdx.y;
	unsigned int bx = blockIdx.x, by = blockIdx.y;
	unsigned int bdx = blockDim.x, bdy = blockDim.y;
	unsigned int col = (bx * bdx) + tx;
	unsigned int row = (by * bdy) + ty;
	unsigned int matrix_pos = row * cols + col;
	bool in_bounds = row < rows && col < cols;

	if(in_bounds) {
		unsigned int dfi_matrix_pos = (row * dfi_matrix_cols) + (dfi_matrix_cols - 1);
		unsigned int max_fj_matrix_pos = (col * max_fj_matrix_cols) + (max_fj_matrix_cols - 1);
		op(output_matrix, input_matrix,
			dfi_matrix[dfi_matrix_pos], max_fj_matrix[max_fj_matrix_pos],
			matrix_pos, cols
		);
	}
}