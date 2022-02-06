#pragma once

#include <cstddef>
#include <vector>
#include <fftw3.h>

class DctCalculator {
public:
    // input and output are width by width matrices, first row, then
    // the second row.
    DctCalculator(size_t width, std::vector<double> *input, std::vector<double> *output);

    void Inverse();

    ~DctCalculator();

private:
    size_t width_;
    std::vector<double> *input_, *output_;
    fftw_plan plan_;
};
