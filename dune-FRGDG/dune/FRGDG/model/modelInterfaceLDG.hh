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

template <typename GV_, unsigned ms0_, unsigned ms1_>
class ModelInterfaceLDGinstat
{
  protected:
    Dune::ParameterTree ptree;

  public:
    using RF = typename GV_::ctype;
    using GV = GV_;

    static constexpr unsigned int LDGlevels = 2;
    static constexpr unsigned int dim = GV_::dimension;
    
    static constexpr unsigned idx = 0;
    static constexpr unsigned ms[2] = { ms0_, ms1_};

    static constexpr unsigned int m = ms[idx];
    
    template<unsigned i_>
      using Range = Dune::FieldVector<RF, ms[i_]>;
    using Range0 = Dune::FieldVector<RF, m>;

    using JacType = Dune::FieldVector<Dune::FieldMatrix<RF,m,m>,dim>;

    RF Lambda;

    RF time;
    RF k;
    RF k2;
    RF k3;
    RF k4;
    RF k5;

    ModelInterfaceLDGinstat(Dune::ParameterTree ptree_)
      : ptree(ptree_), time(0.0)
    {
      Lambda = ptree.get("param.Lambda" , RF(0.7));
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

    static Dune::FieldVector<RF,dim> getEOMIntersection(const Dune::ParameterTree& ptree)
    {
      return Dune::FieldVector<RF,dim>(0.);
    }

    static std::string getName(const Dune::ParameterTree& ptree)
    {
      return "default";
    }
    
    //--------------------------------------------------------------------------------------------------------
    // You may want to overwrite most of these functions in model definitions
 
    Range0 u_max() 
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    Range0 u_min()
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    Range0 has_max()
    {
      return Dune::FieldVector<bool,m>(false);
    }
    Dune::FieldVector<bool,m> has_min()
    {
      return Dune::FieldVector<bool,m>(false);
    }

    // right hand side
    template <typename E, typename X>
      Range0 q(const E &e, const X &x) const
      {
        return Range0(0.0);
      }
    template <typename E, typename X, typename ...R>
      Range0 q(const E &e, const X &x, R... r) const
      {
        return Range0(0.0);
      }

    // initial condition
    template <typename E, typename X>
      Range0 u0(const E &e, const X &x) const
      {
        return Range0(0.);
      }
    
    // jacobian 
    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), JacType>::type
      Jacobian(const E &cell, 
          const X &x,
          const Range<0> &u,
          const Range<1> &g) const
      {
        JacType res(0.);
        return res;
      }

    Range0 apply_R(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    Range0 apply_R_inv(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), void>::type
      max_eigenvalue(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Range<0> &u_s,
          const Range<0> &u_n,
          const Range<1> &g_s,
          const Range<1> &g_n,
          Dune::FieldMatrix<RF, m, dim> &alpha,
          RF &alpha_t) const
      {
      }

    //Flux function
    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), void>::type
      flux(const E &e, const X &x,
          const Range<0> &u, const Range<1> &g,
          Dune::FieldMatrix<RF, m, dim> &F) const
      {
      }

    //Nonconservative Flux function
    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), void>::type
      numericalNonConFlux(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Dune::FieldVector<RF, dim> n_F,
          const Range<0> &u_s, const Range<0> &u_n,
          const Range<1> &g_s, const Range<1> &g_n,
          Range<0> &D_s, Range<0> &D_n) const
      {
        D_s = Range<0>(0.);
        D_n = Range<0>(0.);
      }

    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), void>::type
      numericalDiffFlux(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const Range<0> &u_s, const Range<0> &u_n,
          const Range<1> &g_s, const Range<1> &g_n,
          Dune::FieldMatrix<RF, m, dim> &A, Dune::FieldMatrix<RF, m, dim> &beta, 
          RF &alpha_t) const
      {
        A = Dune::FieldMatrix<RF, m, dim>(0.);
        beta = Dune::FieldMatrix<RF, m, dim>(0.);
      }

    template <typename E, typename X>
      typename std::enable_if<sizeof(E) && (LDGlevels == 2), void>::type
      diffFlux(const E &e, const X &x,
          const Range<0> &u, const Range<1> &g,
          Dune::FieldMatrix<RF, m, dim> &A) const
      {				
        A = Dune::FieldMatrix<RF, m, dim>(0.);
      }	

    // dirichlet boundary condition
    template <typename I, typename X>
      Range0 g(const I &is, const X &x, const Range0 &s) const
      {
        return s;
      }

    // neumann boundary condition for dynamic class
    template<typename I, typename X>
      Range0 j (const I& i, const X& x, const Range0 &s, const Range0 &f) const
      {
        return f;
      }
    template<typename I, typename X>
      Range0 j (const I& i, const X& x, const Range0 &s, const Range<1> &p, const Range0 &f) const
      {
        return f;
      }

    template<typename I, typename X>
      bool b (const I& i, const X& x) const
      {
        return true;
      }
};

template <typename GV_, unsigned ms0_, unsigned ms1_>
class ModelInterfaceLDGstat
{
  protected:
    Dune::ParameterTree ptree;
    
  public:
    using RF = typename GV_::ctype;
    using GV = GV_;

    static constexpr unsigned dim = GV_::dimension;

    static constexpr unsigned idx = 1;
    static constexpr unsigned curIdx = 1;
    static constexpr unsigned depIdx = 0;
    static constexpr unsigned ms[2] = { ms0_, ms1_};

    static constexpr unsigned m = ms[curIdx];
    static constexpr unsigned curm = ms[curIdx];
    static constexpr unsigned depm = ms[depIdx];
    
    template<unsigned i_>
      using Range = Dune::FieldVector<RF, ms[i_]>;
    using curRange = Dune::FieldVector<RF, curm>;
    using depRange = Dune::FieldVector<RF, depm>;

    using JacType = Dune::FieldVector<Dune::FieldMatrix<RF,m,m>,dim>;

    RF Lambda;

    RF time;
    RF k;
    RF k2;
    RF k3;
    RF k4;
    RF k5;

    ModelInterfaceLDGstat(Dune::ParameterTree ptree_)
      : ptree(ptree_), time(0.0)
    {
      Lambda = ptree.get("param.Lambda" , RF(0.7));
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

    static Dune::FieldVector<RF,dim> getEOMIntersection(const Dune::ParameterTree& ptree)
    {
      return Dune::FieldVector<RF,dim>(0.);
    }

    static std::string getName(const Dune::ParameterTree& ptree)
    {
      return "default";
    }
    
    //--------------------------------------------------------------------------------------------------------
    // You may want to overwrite most of these functions in model definitions
 
    curRange u_max() 
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    curRange u_min()
    {
      return Dune::FieldVector<RF,m>(0.);
    }
    curRange has_max()
    {
      return Dune::FieldVector<bool,m>(false);
    }
    Dune::FieldVector<bool,m> has_min()
    {
      return Dune::FieldVector<bool,m>(false);
    }

    // right hand side
    template <typename E, typename X>
      curRange q(const E &e, const X &x) const
      {
        return curRange(0.0);
      }
    template <typename E, typename X, typename ...R>
      curRange q(const E &e, const X &x, R... r) const
      {
        return curRange(0.0);
      }

    // initial condition
    template <typename E, typename X>
      curRange u0(const E &e, const X &x) const
      {
        return curRange(0.);
      }
    
    // jacobian 
    template <typename E, typename X>
      JacType Jacobian(const E &cell, 
          const X &x,
          const depRange &dep,
          const curRange &cur) const
      {
        std::vector<Dune::FieldMatrix<RF,m,m>> res(dim);
        return res;
      }

    curRange apply_R(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    curRange apply_R_inv(const JacType& Jacobian, const unsigned int dir, const Dune::FieldVector<RF,m>& v) const
    {
      Dune::FieldVector<RF,m> v_tr = v;
      return v_tr;
    }

    template <typename E, typename X>
      void max_eigenvalue(const E &inside, const X &x_inside,
          const E &outside, const X &x_outside,
          const depRange &dep_s, const depRange &dep_n,
          const curRange &cur_s, const curRange &cur_n,
          Dune::FieldMatrix<RF, m, dim> &alpha) const
      {
      }

    //Flux function
    template <typename E, typename X>
      void flux(const E &e, const X &x,
          const depRange &dep, const curRange &cur,
          Dune::FieldMatrix<RF, m, dim> &F) const
      {
      }

    // dirichlet boundary condition
    template <typename I, typename X>
      curRange g(const I &is, const X &x, const curRange &s) const
      {
        return s;
      }

    // dirichlet boundary condition
    template <typename I, typename X, typename _Range>
      curRange g(const I &is, const X &x, const _Range &s, const curRange &f) const
      {
        return s;
      }

    // neumann boundary condition for dynamic class
    template<typename I, typename X, typename _Range>
      curRange j(const I& i, const X& x, const _Range &s, const curRange &f) const
      {
        return f;
      }

    template<typename I, typename X>
      bool b (const I& i, const X& x) const
      {
        return true;
      }
};
