//
// Created by Martin Steinegger on 7/23/24.
//

#ifndef FWBW_FWBWALIGNER_H
#define FWBW_FWBWALIGNER_H


#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <string.h>  // for memcpy
#include "SubstitutionMatrix.h"


class FwBwAligner {
public:
    typedef struct {
        uint32_t score1;
        uint32_t score2;
        int32_t dbStartPos1;
        int32_t dbEndPos1;
        int32_t	qStartPos1;
        int32_t qEndPos1;
        int32_t ref_end2;
        float qCov;
        float tCov;
        std::string cigar;
        double evalue;
        int identicalAACnt;
        int32_t cigarLen;
        int word;
    } s_align;

    FwBwAligner(unsigned int maxQueryLen, unsigned int maxTargetLen,
                SubstitutionMatrix & subMat3Di, SubstitutionMatrix & subMat);
    ~FwBwAligner();

    s_align align(std::string & query3Di, std::string & queryAA,
                  std::string & target3Di, std::string & targetAA);

    void computeForwardScoreMatrix(const unsigned char* ss1_num, const unsigned char* ss2_num,
                                   const unsigned char* aa1_num, const unsigned char* aa2_num,
                                   unsigned int queryLen, unsigned int targetLen,
                                   float** mat3di, float** blosum, float T, float ** scoreForward);

    unsigned char* seq2num(const std::string& seq, const unsigned char* aa2num) {
        unsigned char* idx = static_cast<unsigned char*>(malloc(seq.size() * sizeof(unsigned char)));
        for (size_t i = 0; i < seq.size(); ++i) {
            idx[i] = aa2num[seq[i]];
        }
        return idx;
    }

    template <typename T>
    T** malloc_matrix(int dim1, int dim2) {
#define ALIGN_FLOAT 16
#define ICEIL(x_int, fac_int) ((x_int + fac_int - 1) / fac_int) * fac_int
        // Compute mem sizes rounded up to nearest multiple of ALIGN_FLOAT
        size_t size_pointer_array = ICEIL(dim1 * sizeof(T*), ALIGN_FLOAT);
        size_t dim2_padded = ICEIL(dim2 * sizeof(T), ALIGN_FLOAT) / sizeof(T);

        T** matrix = (T**) mem_align(ALIGN_FLOAT, size_pointer_array + dim1 * dim2_padded * sizeof(T));
        if (matrix == NULL)
            return matrix;

        T* ptr = (T*) (matrix + (size_pointer_array / sizeof(T*)));
        for (int i = 0; i < dim1; ++i) {
            matrix[i] = ptr;
            ptr += dim2_padded;
        }
#undef ICEIL
        return matrix;
    }

    inline void* mem_align(size_t boundary, size_t size) {
        void* pointer;
        if (posix_memalign(&pointer, boundary, size) != 0) {
#define MEM_ALIGN_ERROR "mem_align could not allocate memory.\n"
            fwrite(MEM_ALIGN_ERROR, sizeof(MEM_ALIGN_ERROR), 1, stderr);
#undef MEM_ALIGN_ERROR
            exit(3);
        }
        return pointer;
    }

private:

    float** zmForward;
    float** zeForward;
    float** zfForward;
    float** zmBackward;
    float** zeBackward;
    float** zfBackward;
    float** scoreForward;
    float** scoreBackward;
    float** P;

    float** zmaxBlocksMaxForward;
    float** zmaxBlocksMaxBackward;

    float* zmaxForward;
    float* zmaxBackward;
    float *vj;
    float *wj;
    float** mat3di;
    float** blosum;
    const SubstitutionMatrix & subMat3Di;
    const SubstitutionMatrix & subMat;

    void forwardBackwardSaveBlockMaxLocal(float** S, float** z_init, float* vj, float* wj,
                                          float T, float go, float ge,
                                          size_t rows, size_t start, size_t end,
            // output is zm, ze, ze, zmax
                                          float** zm, float** ze, float** zf, float* zmax);

    void rescaleBlocks(float **matrix, float **scale, size_t rows, size_t cols, int length, int blocks);

};


#endif //FWBW_FWBWALIGNER_H
