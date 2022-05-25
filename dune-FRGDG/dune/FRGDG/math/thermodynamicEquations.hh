#pragma once

#include <cmath>

#include <dune/FRGDG/common/utils.hh>

namespace ThermodynamicEquations
{
  static constexpr double epsilon_ = 1e-11;
  static constexpr double pi2 = M_PI*M_PI;

  using utils::powr;
  
  // for convenience, all hyperbolic functions 
  template <typename RF>
    RF Cosh(const RF& x)
    {
      return std::cosh(x);
    }
  template <typename RF>
    RF Sinh(const RF& x)
    {
      return std::sinh(x);
    }
  template <typename RF>
    RF Tanh(const RF& x)
    {
      return std::tanh(x);
    }
  template <typename RF>
    RF Coth(const RF& x)
    {
      return 1./std::tanh(x);
    }
  template <typename RF>
    RF Sech(const RF& x)
    {
      return 1./std::cosh(x);
    }
  template <typename RF>
    RF Csch(const RF& x)
    {
      return 1./std::sinh(x);
    }
  
  // coth(e/2T)
  template <typename RF>
    RF cothS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./std::tanh(e/(T*2.)) : RF(e>0) - RF(e<0);
    }

  // d/de cothS
  template <typename RF>
    RF dcothS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? -1./powr<2,RF>(std::sinh(e/(T*2.)))/(2.*T) : 0.;
    }

  // d^2/de^2 cothS
  template <typename RF>
    RF ddcothS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? - dcothS(e,T)*cothS(e,T)/(2.*T) : 0.;
    }

  // tanh(e/2T)
  template <typename RF>
    RF tanhS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? std::tanh(e/(T*2.)) : RF(e>0) - RF(e<0);
    }

  // d/de tanhS
  template <typename RF>
    RF dtanhS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./powr<2,RF>(std::cosh(e/(T*2.)))/(2.*T) : 0.;
    }

  // d^2/de^2 tanhS
  template <typename RF>
    RF ddtanhS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? - dtanhS(e,T)*tanhS(e,T)/(2.*T) : 0.;
    }

  // sech(e/2T)
  template <typename RF>
    RF sechS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./std::cosh(e/(T*2.)) : 0.;
    }

  // csch(e/2T)
  template <typename RF>
    RF cschS(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./std::sinh(e/(T*2.)) : 0.;
    }
  
  //////////////////////// Distribution Functions nB, nF and their Derivatives

  //Bosonic distribution
  template <typename RF>
  [[deprecated("please use hyperbolic functions instead of the distributions.")]]
    RF nB(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./(std::exp(e/T) - 1.) : RF(e>0.)-1;
    }

  //Derivative d/de of the bosonic distribution
  template <typename RF>
  [[deprecated("please use hyperbolic functions instead of the distributions.")]]
    RF dnB(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? -std::exp(e/T)/powr<2,RF>(std::exp(e/T) - 1.)/T : RF(0.);
    }

  //Derivative d²/de² of the bosonic distribution
  template <typename RF>
  [[deprecated("please use hyperbolic functions instead of the distributions.")]]
    RF ddnB(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? std::exp(e/T)*(1+std::exp(e/T))/powr<3,RF>(std::exp(e/T) - 1.)/powr<2,RF>(T) : RF(0.);
    }

  //Fermionic distribution
  template <typename RF>
  [[deprecated("please use hyperbolic functions instead of the distributions.")]]
    RF nF(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ? 1./(std::exp(e/T) + 1.) : RF(e<0.);
    }

  //Derivative d/de of the fermionic distribution
  template <typename RF>
  [[deprecated("please use hyperbolic functions instead of the distributions.")]]
    RF dnF(const RF& e, const RF& T)
    {
      return (std::abs(T/e)>epsilon_) ?  -std::exp(e/T)/powr<2,RF>(std::exp(e/T) + 1.)/T : RF(0.);
    }
};
