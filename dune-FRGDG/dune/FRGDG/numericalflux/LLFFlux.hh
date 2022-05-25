#pragma once

#include <dune/FRGDG/numericalflux/fluxInterface.hh>

//local Lax-Friedrichs Flux
template <typename MODEL>
class LLFflux : public FluxInterface<MODEL>
{
  public:
    static constexpr unsigned int dim = MODEL::dim;
    static constexpr unsigned int m = MODEL::m;

    using Model = MODEL;
    using RF = typename Model::RF;

    LLFflux(MODEL &model) : model_(model) {}

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
      Dune::FieldMatrix<RF, m, dim> alpha(0.0);
      RF alpha_t(0.0);
      model().max_eigenvalue(inside, x_inside, outside, x_outside, u_s, u_n, alpha, alpha_t);

      for (unsigned i = 0; i < m; i++)
        f[i] = f[i] - 0.5 * alpha[i][std::abs(axis)] * (u_n[i] - u_s[i]);
    
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

//LDG local Lax-Friedrichs Flux
template <typename MODEL>
class LLFfluxLDG : public FluxInterface<MODEL>
{
  public:
    static constexpr unsigned dim = MODEL::dim;
    static constexpr unsigned m0 = MODEL::m0;
    static constexpr unsigned m1 = MODEL::m1;
		
    using Model = MODEL;
    using RF = typename Model::RF; // type for computations
    
    static constexpr unsigned int m = Model::m;

    LLFfluxLDG(MODEL &model) : model_(model) {}

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

			Dune::FieldMatrix<RF, m, dim> A(0.0);
			Dune::FieldMatrix<RF, m, dim> beta(0.0);
      RF alpha_t(0.0);
			
			model().numericalDiffFlux(inside, x_inside, outside, x_outside, u_s, u_n, g_s, g_n, A, beta, alpha_t);
			A.umv(n_F, f);
			
      //max eigenvalue
      Dune::FieldMatrix<RF, m, dim> alpha(0.0);
			model().max_eigenvalue(inside, x_inside, outside, x_outside, u_s, u_n, g_s, g_n, alpha, alpha_t);
			
      //add diffusion
      for (unsigned i = 0; i < m; i++)
        f[i] = f[i] - 0.5 * alpha[i][std::abs(axis)] * (u_n[i] - u_s[i]) - 0.5 * beta[i][std::abs(axis)];
				
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

			Dune::FieldMatrix<RF, m, dim> beta(0.0);

			model().numericalDiffFlux(inside, x_inside, outside, x_outside, u_s, u_n, beta);

			//add LFF
      for (unsigned i = 0; i < m; i++)
        f[i] = f[i] - 0.5 * beta[i][std::abs(axis)];

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
