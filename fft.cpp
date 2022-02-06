
#include "fft.h"
#include <stdexcept>
#include <cmath>

DctCalculator::DctCalculator(size_t width, std::vector<double> *input, std::vector<double> *output)
    : width_(width), input_(input), output_(output) {
    if (!input || !output || input->size() != width * width || output->size() != width * width) {
        throw std::invalid_argument("Invalid matrix size");
    }
    plan_ = fftw_plan_r2r_2d(width_, width_, &input_->front(), &output_->front(), FFTW_REDFT01,
                             FFTW_REDFT01, FFTW_ESTIMATE);
}

DctCalculator::~DctCalculator() {
    fftw_destroy_plan(plan_);
}

void DctCalculator::Inverse() {

    double tmp = 1 / 16.0;
    for (size_t i = 0, k = 0; i != width_; ++i) {
        for (size_t j = 0; j != width_; ++j) {
            double u = i == 0 ? M_SQRT2 : 1;
            double v = j == 0 ? M_SQRT2 : 1;
            (*input_)[k++] *= u * v * tmp;
        }
    }
    fftw_execute(plan_);
}
