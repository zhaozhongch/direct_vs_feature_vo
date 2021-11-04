#include "ExtractHighGradient.cuh"


__global__ void ExtractHighGradientKernel(float* im, float* depth, cv::KeyPoint* kp, float* kp_d, int rows, int cols, int pattern_row, int pattern_col, const int cell, const int pattern_size){
    // //NOTE THIS.......if there are a lot of pixel that index_x is larger than threadIdx.x_max, then they won't be reached
    int thread_idx = threadIdx.x + blockIdx.x * blockDim.x;
    int thread_idy = threadIdx.y + blockIdx.y * blockDim.y;
    //printf("tx, bx %d %d \n", threadIdx.x, blockIdx.x);
    //printf("bm x, dm y %d %d \n", blockDim.x, blockDim.y);
    //int pattern = 3;
    int p = pattern_size;
    int index;
    int cell_length = rows/cell;
    int cell_width = cols/cell;
    float gradient_x, gradient_y, norm;
    float max_norm = 0.0;
    if(thread_idx>pattern_col-1 || thread_idy>pattern_row-1)
        return;
    // if(thread_idy == 0)
    //     printf(" %d %d \n", blockDim.x, blockDim.y);
    for(int i = 0; i < p; i++)
        for(int j = 0; j< p; j++){
            index = i + p * (thread_idx+1) + (p*(thread_idy+1) + j)*cols;
            gradient_x = std::abs(im[index + 1] - im[index - 1]);//right and left
            gradient_y = std::abs(im[index + cols] - im[index - cols]);//down and up
            norm = sqrt(gradient_x*gradient_x + gradient_y*gradient_y);
            //push_to_key points
            if(norm>max_norm){
                kp[thread_idx + thread_idy*pattern_col].pt.y = p*(thread_idy+1) + j;
                kp[thread_idx + thread_idy*pattern_col].pt.x = p*(thread_idx+1) + i;
                kp_d[thread_idx + thread_idy*pattern_col] = depth[index];
                int cell_row = (p*(thread_idy+1) + j)/cell_length;
                int cell_col = (p*(thread_idx+1) + i)/cell_width;
                kp[thread_idx + thread_idy*pattern_col].class_id = cell * cell_row + cell_col;//class_id encodes which cell it belongs to
                max_norm = norm;
            }
        }
    // if(thread_idx == pattern_col-1 && thread_idy == pattern_row-1 && cell == 4 && pattern_size ==7){
    //     printf("final keypoint row and col %f, %f \n", kp[thread_idx + thread_idy*pattern_col].pt.y, kp[thread_idx + thread_idy*pattern_col].pt.x);
    //     printf("it depth and norm is %f, %f \n", max_norm, kp_d[thread_idx + thread_idy*pattern_col]);
    // }
    if(max_norm == 0.0 || kp_d[thread_idx + thread_idy*pattern_col] < 0.001){
        kp[thread_idx + thread_idy*pattern_col].class_id = -1;
    }
};

namespace ORB_SLAM2
{
void ExtractHighGradient(const cv::Mat& im, const cv::Mat& depth, std::vector<cv::KeyPoint>& hg_points, std::vector<float>& hg_depth, const int cell, const int pattern_size){
    // cv::Mat test;
    // std::string add("/home/zhaozhong/dataset/tum/rgbd_dataset_freiburg1_teddy/rgb/1305032180.416768.png");//image read will cost about 6~7 ms
    // test = cv::imread(add, 0);//0 for gray and 1 for rgb

    int rows = im.rows;
    int cols = im.cols;
    int p = pattern_size;
    int row_pattern_num = (rows - 2*p) / p;
    int col_pattern_num = (cols - 2*p) / p;
    int keypoint_num = row_pattern_num * col_pattern_num;
    //std::cout<<"col pattern num and row pattern num "<<col_pattern_num<<", "<<row_pattern_num<<std::endl;
    //std::cout<<"keypoint num is "<<keypoint_num<<std::endl;
    int image_size = sizeof(float) * rows * cols;
    int depth_image_size = sizeof(float) * rows * cols;
    int keypoint_size = sizeof(cv::KeyPoint) * keypoint_num;
    int depth_size = sizeof(float) * keypoint_num;
    
    //printf("size of keypoint %ld \n", sizeof(cv::KeyPoint));

    //double extract_time = (double)cv::getTickCount();
    //1 allocate memory in the host(CPU)
    cv::KeyPoint* host_output = (cv::KeyPoint*)malloc(keypoint_size);
    float* host_output_d = (float*)malloc(depth_size);

    //2 assign value to the data (In reality, you get your data from a thousand ways)
    //host_input = test.data;

    //3 allocate memory in the device (GPU)
    float* device_input;
    float* device_input_d;
    cv::KeyPoint* device_output;
    float* device_output_d;
    
    //double alloc_time = (double)cv::getTickCount();
    cudaMalloc((void**)&device_output, keypoint_size);
    cudaMalloc((void**)&device_input, image_size);
    cudaMalloc((void**)&device_output_d, depth_size);
    cudaMalloc((void**)&device_input_d, depth_image_size);
    
    //4 copy the data from the host to the device
    //1 byte element size for uchar, 4 byte for float
    cudaMemcpy(device_input, im.ptr<float>(0), image_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_input_d, depth.ptr<float>(0), depth_image_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_output, host_output, keypoint_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_output_d, host_output_d, depth_size, cudaMemcpyHostToDevice);

    //5 assign block number and thread number in the block
    // On my lenovo computer
    // Maximum number of threads per multiprocessor:  2048
    // Maximum number of threads per block:           1024
    // 6 multi processor
    dim3 block_number = dim3(3, 3);
    dim3 thread_number = dim3(col_pattern_num/3 + 1, row_pattern_num/3 + 1);//(blockDim.x, blockDim.y), x for col y ,for row,
    //std::cout<<"thread number in one block "<<col_pattern_num/3 + 1<<", "<<row_pattern_num/3 + 1<<std::endl;

    //6 call function in the device. During this step, results should be ready.
    ExtractHighGradientKernel<<<block_number, thread_number>>>(device_input, device_input_d, device_output, device_output_d, rows, cols, row_pattern_num, col_pattern_num, cell, p);//in cpp file '<<<' doesn't make sense and will lead to error

    //7 copy memory from device to host
    cudaMemcpy(host_output, device_output, keypoint_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_output_d, device_output_d, depth_size, cudaMemcpyDeviceToHost);

    for(int i = 0; i<keypoint_num; i++){
        if(host_output[i].class_id != -1){
            hg_points.push_back(host_output[i]);
            hg_depth.push_back(host_output_d[i]);
        }
    }
    //std::cout<<"final point size is "<<hg_points.size()<<std::endl;

    // extract_time = ((double)cv::getTickCount() - extract_time)/cv::getTickFrequency();
    // std::cout<<"extract time is "<<extract_time<<std::endl;

    //std::vector<cv::KeyPoint> hg_points(host_output, host_output + keypoint_num)

    // if(p==7 && cell==4){
    //     cv::Mat im_with_keypoints;

    //     cv::Mat show_case = im.clone();
    //     cv::Mat show_case2;
    //     show_case.convertTo(show_case2, CV_8U);
    //     cv::drawKeypoints(show_case2, hg_points, im_with_keypoints);//draw image here to show poitns may lead to error or not. Just use this to visualize for a moment

    //     cv::namedWindow("show key points", cv::WINDOW_AUTOSIZE);
    //     cv::imshow("show key points", im_with_keypoints);

    //     cv::waitKey(1);
    // }

    //8 free the momory in device and the host
    //free(test.data);//host_input
    free(host_output);//free host_output won't influence the value in the hg_points
    free(host_output_d);
    cudaFree(device_input_d);
    cudaFree(device_output_d);
    cudaFree(device_input);
    cudaFree(device_output);
}
}