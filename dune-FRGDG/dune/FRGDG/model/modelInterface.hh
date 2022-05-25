#pragma once

// this module
#include <dune/FRGDG/common/utils.hh>

#include <dune/FRGDG/solver/ovlpistlsolverbackend_extension.hh>
#include <dune/FRGDG/integrationscheme/LDGScheme.hh>
#include <dune/FRGDG/integrationscheme/DGScheme.hh>

// dune includes
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/finiteelementmap/qkdg.hh>

// std
#include <string>

template <typename GV_, unsigned int m_ = 1>
class ModelInterfaceConEq
{
  protected:
    Dune::ParameterTree ptree;

  public:
    using RF = typename GV_::ctype;
    using GV = GV_;
    static constexpr unsigned dim = GV_::dimension;
    static constexpr unsigned m = m_;
    using Range = Dune::FieldVector<RF, m>;

    using JacType = Dune::FieldVector<Dune::FieldMatrix<RF,m,m>,dim>;

    RF time;
    RF k;
    RF k2;
    RF k3;
    RF k4;
    RF k5;
    RF Lambda;

    ModelInterfaceConEq(Dune::ParameterTree ptree_)
      : ptree(ptree_), time(0.0) 
    {
      Lambda = ptree.get("param.Lambda" , RF(0.7));
    }

    static Dune::FieldVector<RF,dim> getEOMIntersection(const Dune::ParameterTree& ptree)
    {
      return Dune::FieldVector<RF,dim>(0.);
    }

    static std::string getName(const Dune::ParameterTree& ptree)
    {
      return "default";
    }

    //--------------------------------------------------------------------------------------------------------
    // These functions are standard utility and access functions
    
    void setTime(RF t)
    {
      time = t;
      k = std::exp(-time) * Lambda;
      k2 = k*k;
      k3 = k2*k;
      k4 = k2*k2;
      k5 = k3*k2;
    }

    const RF& getTime() const
    {
      return time;
    }

    const RF& getk2() const
    {
      return k2;
    }
    
    //--------------------------------------------------------------------------------------------------------
    // You may want to overwrite most of these functions in model definitions
 
    Dune::FieldVector<RF,m> u_max() 
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    Dune::FieldVector<RF,m> u_min()
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    Dune::FieldVector<bool,m> has_max()
    {
      return Dune::FieldVector<bool,m>(false);
    }
    Dune::FieldVector<bool,m> has_min()
    {
      return Dune::FieldVector<bool,m>(false);
    }

    // right hand side
    template <typename E, typename X>
      Range q(const E &e, const X &x) const
      {
        return Range(0.0);
      }
    template <typename E, typename X, typename ...R>
      Range q(const E &e, const X &x, R... r) const
      {
        return Range(0.0);
      }

    // initial condition
    template <typename E, typename X>
      Range u0(const E &e, const X &x) const
      {
        return Range(0.);
      }
    
    // jacobian 
    template <typename E, typename X>
      JacType Jacobian(const E &cell, 
          const X &x,
          const Dune::FieldVector<RF, m> &u) const
      {
        JacType res(0.);
        return res;
      }

    Dune::FieldVector<RF,m> apply_R(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    Dune::FieldVector<RF,m> apply_R_inv(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    template <typename E, typename X>
      void max_eigenvalue(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Dune::FieldVector<RF, m> &u_s,
          const Dune::FieldVector<RF, m> &u_n,
          Dune::FieldMatrix<RF, m, dim> &alpha,
          RF &alpha_t) const
      {
      }

    //Flux function
    template <typename E, typename X>
      void flux(const E &e, const X &x,
          const Dune::FieldVector<RF, m> &u,
          Dune::FieldMatrix<RF, m, dim> &F) const
      {
      }

    // dirichlet boundary condition
    template <typename I, typename X>
      Range g(const I &is, const X &x, const Range &s) const
      {
        return Range(s);
      }

    // dirichlet boundary condition
    template <typename I, typename X>
      Range g(const I &is, const X &x, const Range &s, const Dune::FieldVector<RF, m> &f) const
      {
        return Range(s);
      }

    // neumann boundary condition
    template<typename I, typename X>
      Range j (const I& i, const X& x, const Range &s, const Dune::FieldVector<RF, m> &f) const
      {
        return Range(f);
      }

    template<typename I, typename X>
      bool b (const I& i, const X& x) const
      {
        return true;
      }
};
