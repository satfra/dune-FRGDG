#pragma once

#include <dune/FRGDG/numericalflux/fluxInterface.hh>

//Central Flux
template <typename MODEL>
class CentralFlux : public FluxInterface<MODEL>
{
  public:
    static constexpr unsigned int dim = MODEL::dim;
    static constexpr unsigned int m = MODEL::m;

    using Model = MODEL;
    using RF = typename Model::RF; // type for computations

    CentralFlux(MODEL &model) : model_(model) {}
    
    // For a Conservation equation
    template <typename E, typename X>
    void numericalFlux(const E &inside, const X &x_inside,
        const E &outside, const X &x_outside,
        const Dune::FieldVector<RF, dim> n_F,
        const Dune::FieldVector<RF, m> &u_s,
        const Dune::FieldVector<RF, m> &u_n,
        const int axis,
        Dune::FieldVector<RF, m> &f, RF &amax) const
    {
      f = 0.0;

      //evaluate flux
      model().flux(inside, x_inside, u_s, Fs);
      model().flux(outside, x_outside, u_n, Fn);

      //Fs*n_F + Fn*n_F
      Fs.umv(n_F, f);
      Fn.umv(n_F, f);
      f *= 0.5;

      //max eigenvalue
      Dune::FieldMatrix<RF, dim, m> alpha(0.0);
      RF alpha_t(0.0);
      model().max_eigenvalue(inside, x_inside, outside, x_outside, u_s, u_n, alpha, alpha_t);

      amax = std::max(alpha_t, amax);
    }

    MODEL &model() const
    {
      return model_;
    }

  private:
    mutable Dune::FieldMatrix<RF, m, dim> Fs;
    mutable Dune::FieldMatrix<RF, m, dim> Fn;
    MODEL &model_;
};

//LDG Central Flux
template <typename MODEL>
class CentralFluxLDG : public FluxInterface<MODEL>
{
  public:
    static constexpr unsigned int dim = MODEL::dim;
    static constexpr unsigned int m = MODEL::m;

    static constexpr unsigned int m0 = MODEL::m0;
    static constexpr unsigned int m1 = MODEL::m1;

    using Model = MODEL;
    using RF = typename Model::RF; // type for computations

    CentralFluxLDG(MODEL &model) : model_(model) {}

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

      //evaluate flux
      model().flux(inside, x_inside, u_s, g_s, Fs);
      model().flux(outside, x_outside, u_n, g_n, Fn);

      //Fs*n_F + Fn*n_F
      Fs.umv(n_F, f);
      Fn.umv(n_F, f);
      f *= 0.5;

      //max eigenvalue
      Dune::FieldMatrix<RF, dim, m> alpha(0.0);
      RF alpha_t(0.0);
      model().max_eigenvalue(inside, x_inside, outside, x_outside, u_s, u_n, g_s, g_n, alpha, alpha_t);
      
      amax = std::max(alpha_t, amax);
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

      //evaluate flux
      model().flux(inside, x_inside, u_s, Fs);
      model().flux(outside, x_outside, u_n, Fn);

      //Fs*n_F + Fn*n_F
      Fs.umv(n_F, f);
      Fn.umv(n_F, f);
      f *= 0.5;
    }

    MODEL &model() const
    {
      return model_;
    }

  private:
    mutable Dune::FieldMatrix<RF, m, dim> Fs;
    mutable Dune::FieldMatrix<RF, m, dim> Fn;
    MODEL &model_;
};
