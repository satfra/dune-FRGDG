#pragma once

#include <dune/FRGDG/numericalflux/fluxInterface.hh>

//LDG Upwind Flux
template <typename MODEL>
class DownwindFluxLDG : public FluxInterface<MODEL>
{
  public:
    static constexpr unsigned int dim = MODEL::dim;
    static constexpr unsigned int m = MODEL::m;
    static constexpr unsigned int m0 = MODEL::m0;
    static constexpr unsigned int m1 = MODEL::m1;

    using Model = MODEL;
    using RF = typename Model::RF; // type for computations

    DownwindFluxLDG(MODEL &model) : model_(model) {}

    // for LDG with one stationary equation
    template <typename E, typename X, typename Range0, typename Range1>
    void numericalFlux(const E &inside, const X &x_inside,
        const E &outside, const X &x_outside,
        const Dune::FieldVector<RF, dim> &n_F,
        const Range0 &u_s, const Range0 &u_n,
        const Range1 &g_s, const Range1 &g_n,
        const int &axis,
        Dune::FieldVector<RF, m> &f, RF &amax) const
    {
      f = 0.0;

      auto pos_s = inside.geometry().center();
      auto pos_n = outside.geometry().center();
      if(pos_s[0] < pos_n[0])
        model().flux(outside, x_outside, u_n, g_n, F);
      else
        model().flux(inside, x_inside, u_s, g_s, F);

      F.umv(n_F, f);
    }

    // for stationary equations
    template <typename E, typename X, typename Range>
    void numericalFlux(const E &inside, const X &x_inside,
        const E &outside, const X &x_outside,
        const Dune::FieldVector<RF, dim> &n_F,
        const Range &u_s, const Range &u_n,
        const int &axis,
        Dune::FieldVector<RF, m> &f) const
    {
      f = 0.0;

      auto pos_s = inside.geometry().center();
      auto pos_n = outside.geometry().center();
      if(pos_s[0] < pos_n[0])
        model().flux(outside, x_outside, u_n, F);
      else
        model().flux(inside, x_inside, u_s, F);

      F.umv(n_F, f);
    }

    MODEL &model() const
    {
      return model_;
    }

  private:
    mutable Dune::FieldMatrix<RF, m, dim> F;
    MODEL &model_;
};
