#pragma once

#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>

namespace Dune {
  namespace PDELab {
    /** @brief This solver will solve the system diag(A)z = r for z exactly on an overlapping grid.
     *
     * @tparam GFS The Type of the GridFunctionSpace.
     */
    template<class GFS>
    class ISTLBackend_OVLP_Diagonal
          : public OVLPScalarProductImplementation<GFS>, public LinearResultStorage
    {
    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
      */
      explicit ISTLBackend_OVLP_Diagonal (const GFS& gfs_)
        : OVLPScalarProductImplementation<GFS>(gfs_), gfs(gfs_)
      {}

      explicit ISTLBackend_OVLP_Diagonal (const ISTLBackend_OVLP_Diagonal& other_)
        : gfs(other_.gfs)
      {}

      /*! \brief solve the given linear system
       
        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename Dune::template FieldTraits<typename W::ElementType >::real_type reduction)
      {
        using Backend::Native;
        using Backend::native;
        Dune::SeqJac<
          Native<M>,
          Native<V>,
          Native<W>
          > jac(native(A),1,1.0);
        jac.pre(native(z),native(r));
        jac.apply(native(z),native(r));
        jac.post(native(z));
        if (gfs.gridView().comm().size()>1)
        {
          CopyDataHandle<GFS,V> copydh(gfs,z);
          gfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = static_cast<double>(reduction);
        res.conv_rate  = static_cast<double>(reduction); // pow(reduction,1.0/1)
      }

    private:
      const GFS& gfs;
    };
  }
}

