//
// Created by Martin Steinegger on 7/23/24.
//

#include "FwBwAligner.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip> // for std::setprecision
#include "simd.h"

FwBwAligner::FwBwAligner(size_t queryLen, size_t targetLen, size_t length, size_t blocks,
                         SubstitutionMatrix &subMat3Di ,SubstitutionMatrix &subMat):
                         subMat3Di(subMat3Di), subMat(subMat){


    zmForward = malloc_matrix<float>(queryLen, targetLen);
    //zeForward = malloc_matrix<float>(queryLen, targetLen);
    //zfForward = malloc_matrix<float>(queryLen, targetLen);
    zmBackward = malloc_matrix<float>(queryLen, targetLen);
    //zeBackward = malloc_matrix<float>(queryLen, targetLen);
    //zfBackward = malloc_matrix<float>(queryLen, targetLen);
    scoreForward = malloc_matrix<float>(queryLen, targetLen);
    scoreBackward = malloc_matrix<float>(queryLen, targetLen);
    P = malloc_matrix<float>(queryLen, targetLen);
    btMatrix = static_cast<uint8_t*>(mem_align(16, queryLen * targetLen * sizeof(uint8_t)));

    // Define block matrices
    // int length = targetLen / 4;
    vj = static_cast<float*>(mem_align(16, length * sizeof(float)));
    wj = static_cast<float*>(mem_align(16, length * sizeof(float)));

    // Matrices used outside forwardBackwardSaveBlockMaxLocal, so shape is blocks, queryLen
    zmaxBlocksMaxForward = malloc_matrix<float>(blocks, queryLen);
    zmaxBlocksMaxBackward = malloc_matrix<float>(blocks, queryLen);

    // Matrices used inside forwardBackwardSaveBlockMaxLocal, so shape is queryLen + 1
    zmaxForward = static_cast<float*>(malloc((queryLen + 1) * sizeof(float)));
    memset(zmaxForward, 0, (queryLen + 1) * sizeof(float)); 
    zmaxBackward = static_cast<float*>(malloc((queryLen + 1) * sizeof(float)));
    memset(zmaxBackward, 0, (queryLen + 1) * sizeof(float));

    zeBlock = static_cast<float*>(malloc((length + 1) * sizeof(float)));
    zfBlock = static_cast<float*>(malloc((length + 1) * sizeof(float)));
    zmBlock = malloc_matrix<float>(queryLen + 1, length + 1);


    mat3di = malloc_matrix<float>(21, 21);
    blosum = malloc_matrix<float>(21, 21);
    
    for (size_t i = 0; i < subMat.alphabetSize; ++i) {
        for (size_t j = 0; j < subMat.alphabetSize; ++j) {
            mat3di[i][j] = static_cast<float>(subMat3Di.subMatrix[i][j]) * 2;
            blosum[i][j] = static_cast<float>(subMat.subMatrix[i][j]) * 2;
        }
    }
}


void FwBwAligner::computeForwardScoreMatrix(const unsigned char* ss1_num, const unsigned char* ss2_num,
                                                const unsigned char* aa1_num, const unsigned char* aa2_num,
                                                unsigned int queryLen, unsigned int targetLen,
                                                float** mat3di, float** blosum, float T, float ** scoreForward) {
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            scoreForward[i][j] = mat3di[ss1_num[i]][ss2_num[j]] + blosum[aa1_num[i]][aa2_num[j]];
            scoreForward[i][j] = exp(scoreForward[i][j] / T);
        }
    }
}

FwBwAligner::~FwBwAligner(){
    free(scoreForward);
    free(scoreBackward);
    free(zmForward);
    //free(zeForward);
    //free(zfForward);
    free(zmBackward);
    //free(zeBackward);
    //free(zfBackward);
    free(P);
    free(btMatrix);
    free(zmaxBlocksMaxForward);
    free(zmaxBlocksMaxBackward);
    free(zmaxForward);
    free(zmaxBackward);
    free(vj);
    free(wj);
    free(blosum);
    free(mat3di);
    free(zeBlock);
    free(zfBlock);
    free(zmBlock);
}


FwBwAligner::s_align FwBwAligner::align(std::string & query3Di, std::string & queryAA,
                                        std::string & target3Di, std::string & targetAA,
                                        size_t length, size_t blocks){
    const unsigned int queryLen = query3Di.size();
    const unsigned int targetLen = target3Di.size();

    unsigned char* ss1_num = seq2num(query3Di, subMat3Di.aa2num);
    unsigned char* ss2_num = seq2num(target3Di, subMat3Di.aa2num);
    unsigned char* aa1_num = seq2num(queryAA, subMat.aa2num);
    unsigned char* aa2_num = seq2num(targetAA, subMat.aa2num);

    const float T = 10;
    const float go = -3.5;
    const float ge = -0.3;
    computeForwardScoreMatrix(ss1_num, ss2_num, aa1_num, aa2_num,
                              queryLen, targetLen, mat3di, blosum, T, scoreForward);
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            scoreBackward[i][j] = scoreForward[queryLen - 1 - i][targetLen - 1 - j];
        }
    }



    /*for (size_t i = 0; i < length; ++i) {
        vj[i] = exp(((length - 1) * ge + go - i * ge) / T);
        wj[i] = exp(((length - 1) * ge - i * ge) / T);
    }*/

    /*for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            zmForward[i][j] = -DBL_MAX;
            zeForward[i][j] = -DBL_MAX;
            zfForward[i][j] = -DBL_MAX;
            zmBackward[i][j] = -DBL_MAX;
            zeBackward[i][j] = -DBL_MAX;
            zfBackward[i][j] = -DBL_MAX;
        }
    }*/
   for (size_t j = 0; j < targetLen; ++j) {
            zmForward[0][j] = -DBL_MAX;
            //zeForward[0][j] = -DBL_MAX;
            //zfForward[0][j] = -DBL_MAX;
            zmBackward[0][j] = -DBL_MAX;
            //zeBackward[0][j] = -DBL_MAX;
            //zfBackward[0][j] = -DBL_MAX;
    }
    for(size_t i = 0; i < queryLen; ++i){
        zmForward[i][0] = -DBL_MAX;
        //zeForward[i][0] = -DBL_MAX;
        //zfForward[i][0] = -DBL_MAX;
        zmBackward[i][0] = -DBL_MAX;
        //zeBackward[i][0] = -DBL_MAX;
        //zfBackward[i][0] = -DBL_MAX;
    }




    //initialize zInit
    float* zInitForward[3];
    zInitForward[0] = new float[queryLen];
    zInitForward[1] = new float[queryLen];
    zInitForward[2] = new float[queryLen];


    for (unsigned int i=0; i < queryLen; ++i){
        zInitForward[0][i] = -DBL_MAX;
        zInitForward[1][i] = -DBL_MAX;
        zInitForward[2][i] = -DBL_MAX;
    }
    //float** zmaxBlocksMaxForward = malloc_matrix<float>(blocks+1, queryLen + 1);
    // memcpy(zmaxBlocksMaxForward[0], zmaxForward, (queryLen + 1) * sizeof(float));
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        // std::cout << "Block: " << b << " start " << start << " end " << end << std::endl;
        // number of cols to memcpy in forwardBackwardSaveBlockMaxLocal
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;
        forwardBackwardSaveBlockMaxLocal(scoreForward, zInitForward, vj, wj, T, go, ge, queryLen, start, end, memcpy_cols,
                                         zmForward, zmaxForward, zmBlock ,zeBlock, zfBlock);
        
        memcpy(zmaxBlocksMaxForward[b], zmaxForward, queryLen * sizeof(float));
    }
    


    ///////////////////////////////////Backward////////////////////////////////////////
    
    float* zInitBackward[3];
    zInitBackward[0] = new float[queryLen];
    zInitBackward[1] = new float[queryLen];
    zInitBackward[2] = new float[queryLen];

    for (unsigned int i=0; i < queryLen; ++i){
        zInitBackward[0][i] = -DBL_MAX;
        zInitBackward[1][i] = -DBL_MAX;
        zInitBackward[2][i] = -DBL_MAX;
    }


    // float** zmaxBlocksMaxBackward = malloc_matrix<float>(blocks+1, queryLen + 1);
    // memcpy(zmaxBlocksMaxBackward[0], zmaxBackward, (queryLen + 1) * sizeof(float));
    // std::cout << "Backward start" << std::endl;
    for (size_t b = 0; b < blocks; ++b) {
        // std::cout << "Block " << b << std::endl;
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;
        
        forwardBackwardSaveBlockMaxLocal(scoreBackward, zInitBackward, vj, wj, T, go, ge, queryLen, start, end, memcpy_cols,
                                         zmBackward, zmaxBackward, zmBlock, zeBlock, zfBlock);
        memcpy(zmaxBlocksMaxBackward[b], zmaxBackward, queryLen * sizeof(float));
    }

    ///////////////////////////////////Rescale////////////////////////////////////////
    // Rescale the values by the maximum in the log space for each block
    // This turns the matrix into log space
    
    rescaleBlocks(zmForward, zmaxBlocksMaxForward, queryLen, length, blocks, targetLen);
    rescaleBlocks(zmBackward, zmaxBlocksMaxBackward, queryLen, length, blocks, targetLen);


    float max_zm = -std::numeric_limits<float>::max();
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            max_zm = std::max(max_zm, zmForward[i][j]);
        }
    }

    float sum_exp= 0.0;
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            sum_exp += exp(zmForward[i][j] - max_zm);
        }
    }
    float logsumexp_zm = max_zm + log(sum_exp);



    /*std::ofstream outfile;
    outfile.open("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < queryLen; ++i) {
            for (size_t j = 0; j < targetLen; ++j) {
                outfile << zmForward[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    } else {
        std::cout << "Unable to open file";
    }*/
    
    // compute posterior probabilities
    float max_p = 0.0;
    size_t max_i;
    size_t max_j;
    float maxP = -std::numeric_limits<float>::max();
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            // Debug(Debug::INFO) << zmForward[i][j] << '\t' << zmBackward[queryLen - 1 - i][targetLen - 1 - j] << '\n';
            P[i][j] = exp(
                zmForward[i][j]
                + zmBackward[queryLen - 1 - i][targetLen - 1 - j]
                - log(scoreForward[i][j]) // FIXME scoreForward is already exp(S/T)
                - logsumexp_zm
            );
            maxP = std::max(maxP, P[i][j]);
            // Debug(Debug::INFO) << P[i][j] << '\t';
        }
        // Debug(Debug::INFO) << '\n';
    }
    std::cout << "maxP\t" << maxP << "\n";
    /*outfile.open("/home/lasse/Desktop/Projects/FB_martin/P_mat.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < queryLen; ++i) {
            for (size_t j = 0; j < targetLen; ++j) {
                outfile << P[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    } else {
        std::cout << "Unable to open file";
    }
    std::cout << max_p << std::endl;*/

    // traceback 
    s_align result;
    /*result.cigar = "";
    result.cigar.reserve(queryLen + targetLen);
    result.qEndPos1 = max_i;
    result.dbEndPos1 = max_j;
    float d;
    float l;
    float u;
    size_t index = 0;
    while (max_i > 0 && max_j > 0) {
        d = P[max_i - 1][max_j - 1];
        l = P[max_i][max_j - 1];
        u = P[max_i - 1][max_j];
        // std::cout << std::fixed << std::setprecision(8) << max_i << '\t' << max_j << '\t' << d << '\t' << l << '\t' << u << '\n';
        if (d > l && d > u) {
            max_i--;
            max_j--;
            result.cigar.push_back('M');
        } else if (l > d && l > u) {
            max_j--;
            result.cigar.push_back('I');
        } else {
            max_i--;
            result.cigar.push_back('D');
        }
    }
    result.qStartPos1 = max_i;
    result.dbStartPos1 = max_j;
    result.cigarLen = result.cigar.length();
    std::reverse(result.cigar.begin(), result.cigar.end());

    free(ss1_num);
    free(ss2_num);
    free(aa1_num);
    free(aa2_num);
    delete[] zInitForward[0];
    delete[] zInitForward[1];
    delete[] zInitForward[2];
    delete[] zInitBackward[0];
    delete[] zInitBackward[1];
    delete[] zInitBackward[2];
    */
    return result;
}

void FwBwAligner::forwardBackwardSaveBlockMaxLocal(float** S, float** z_init,
                                                   float* vj, float* wj,
                                                   float T, float go, float ge,
                                                   size_t rows, size_t start, size_t end, size_t memcpy_cols,
                                                   float** zm, float* zmax, float** zmBlock, float* zeBlock, float* zfBlock) {
    float exp_go = exp(go / T);
    float exp_ge = exp(ge / T);
    
    

    //std::cout << "test" << std::endl;
    memset(zeBlock, 0, (end - start + 1) * sizeof(float)); 
    memset(zfBlock, 0, (end - start + 1) * sizeof(float)); 
 
    std::vector<float> ze_first(rows+1, 0);
    std::vector<float> zf_first(rows+1, 0);
    
    //Init blocks
    memset(zmBlock[0], 0, (end - start + 1) * sizeof(float));

    for (size_t i = 0; i < rows; ++i) {
        zmBlock[i+1][0] = z_init[0][i];
        ze_first[i+1] = z_init[1][i];
        zf_first[i+1] = z_init[2][i];
    }


    size_t cols = memcpy_cols;
    float* exp_ge_arr = static_cast<float*>(mem_align(16, cols * sizeof(float)));
    //for (size_t i = 0; i < cols; ++i) {
    //    exp_ge_arr[i] = exp((i * ge + ge) / T);
    //}

    float current_max = 0;
    for (size_t i = 1; i <= rows; ++i) {
        if (i != 1) {
            zmBlock[i - 1][0] = exp(zmBlock[i - 1][0]);
            ze_first[i - 1] = exp(ze_first[i - 1]);
            zf_first[i - 1] = exp(zf_first[i - 1]);
            // Debug(Debug::INFO) << zmBlock[i - 1][0] << '\t';
        }
        const float expMax = exp(-current_max);
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            if (j == 1){
                float tmp = (zmBlock[i - 1][j - 1] + ze_first[i-1] + zf_first[i - 1] + expMax);
                zmBlock[i][j] = tmp * S[i - 1][start + j - 1];
            }
            else{
                float tmp = (zmBlock[i - 1][j - 1] + zeBlock[j - 1] + zfBlock[j - 1] + expMax);
                zmBlock[i][j] = tmp * S[i - 1][start + j - 1];
            }
        }
        

        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            if (j == 1) {
                zeBlock[j] = exp(zmBlock[i][j - 1]) * exp_go + exp(ze_first[i]) * exp_ge;
            } else {
                zeBlock[j] = zmBlock[i][j - 1] * exp_go + zeBlock[j - 1] * exp_ge;
            }

        }
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
                zfBlock[j] = zmBlock[i - 1][j] * exp_go + zfBlock[j] * exp_ge;
        }

        float z_temp = *std::max_element(zmBlock[i] + 1, zmBlock[i] + cols + 1);
        zmax[i-1] = log(z_temp);
        //zmax[i-1] = z_temp;
        current_max += zmax[i-1];
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            zmBlock[i][j] /= z_temp;
            zeBlock[j] /= z_temp;
            zfBlock[j] /= z_temp;
        }
        
        
        zmBlock[i][0] -= zmax[i-1];
        ze_first[i] -= zmax[i-1];
        zf_first[i] -= current_max;
        // zfBlock[i][0] -= zmax[i-1];
        if (i < rows) {
            zmBlock[i+1][0] -= current_max;
            ze_first[i+1] -= current_max;
            z_init[0][i-1] = log(zmBlock[i][cols]) + current_max;
            z_init[1][i-1] = log(zeBlock[cols]) + current_max;
            z_init[2][i-1] = log(zfBlock[cols]) + current_max;

        }
        
        // for (size_t j = i; j <= rows; ++j) {
        //     zmBlock[j][0] -= zmax[i-1];
        //     zeBlock[j][0] -= zmax[i-1];
        //     zfBlock[j][0] -= zmax[i-1];
        // }

        //free(zm_exp);
    }
    // Debug(Debug::INFO) << '\n';

    //Calculate the cumulative sum of zmax[1:]
    std::vector<float> rescale(rows);
    // std::partial_sum(zmax + 1, zmax + rows + 1, rescale.begin());
    std::partial_sum(zmax, zmax + rows, rescale.begin());

    //Fixme
    // 
    /*for (size_t i = 0; i < rows; ++i) {
        z_init[0][i] = log(zmBlock[i + 1][memcpy_cols]) + rescale[i];
        z_init[1][i] = log(zeBlock[i + 1][memcpy_cols]) + rescale[i];
        z_init[2][i] = log(zfBlock[i + 1][memcpy_cols]) + rescale[i];
    }*/

    for (size_t i = 0; i < rows; ++i) {
        memcpy(zm[i] + start, zmBlock[i+1]+1, memcpy_cols * sizeof(float));
        //memcpy(ze[i] + start, zeBlock[i+1]+1, memcpy_cols * sizeof(float));
        //memcpy(zf[i] + start, zfBlock[i+1]+1, memcpy_cols * sizeof(float));
    }

    //free(zmBlock);
    //free(zeBlock);
    //free(zfBlock);
    //free(exp_ge_arr);
}

void FwBwAligner::rescaleBlocks(float **matrix, float **scale, size_t rows, size_t length, size_t blocks, size_t targetLen){
    // Function to rescale the values by the maximum in the log space for each block
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = std::min((b + 1) * length, targetLen);
        // size_t end = (b + 1) * length;
        std::vector<float> cumsum(rows);
        std::partial_sum(scale[b], scale[b] + rows, cumsum.begin());
        // DEBUG:: print cumsum vector for each block
        // std::cout << "block " << b << " cumsum: ";
        // for (size_t i = 0; i < rows; ++i) {
        //     std::cout << cumsum[i] << " ";
        // }
        // std::cout << std::endl;

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = start; j < end; ++j) {
                // matrix[i][j] = log(matrix[i][j]) ;// + cumsum[i-1];
                matrix[i][j] = log(matrix[i][j]) + cumsum[i];
            }
        }
    }
}
