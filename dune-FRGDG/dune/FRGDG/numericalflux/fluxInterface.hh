#pragma once

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

template <typename MODEL_>
class FluxInterface
{
  public:
    using MODEL = MODEL_;
    static constexpr unsigned int dim = MODEL::dim;
    static constexpr unsigned int m = MODEL::m;

    FluxInterface() {};

    using Model = MODEL;
    using RF = typename Model::RF; // type for computations
};
