// COPYRIGHT

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#define protected public  // heh heh
#include <complex>
#include <mtl/mtl.h>
#include <mtl/matrix.h>
#include <mtl/blais.h>
#include <mtl/dense1D.h>
#include "mtl_specializations.h"

////
// Implemententation of  matrix specializations using atlas/blais functionality
////
#ifdef HAVE_BLAS
#include "BlasHeaders.h"

/* minimal M*N*K to use cblas routines */
#define N_CBLAS 2000

namespace mtl {
  
  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * ONE);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(ONE * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * ONE);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const double *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const double *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      double *m3_d=m3.data();
      double fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        double *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const double *m1p=m1_d+r*m1_rs;
          const double *m2p=m2_d+c*m2_cs;
          double tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const double ONE=double(1.0);
      const double alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=m1.ncols();
      size_t m1_cs=1;
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.ncols(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * ONE);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=1.0*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(ONE * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=m2.ncols();
      size_t m2_cs=1;
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.ncols(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*1.0;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * ONE);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=m3.ncols();
      size_t m3_cs=1;
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.ncols());
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    if (nr*nc*k < N_CBLAS) {
      size_t m1_rs=1;
      size_t m1_cs=m1.nrows();
      const float *m1_d=m1.get_twod().twod.data();
      size_t m2_rs=1;
      size_t m2_cs=m2.nrows();
      const float *m2_d=m2.get_twod().twod.data();
      size_t m3_rs=1;
      size_t m3_cs=m3.nrows();
      float *m3_d=m3.data();
      float fac=m1.get_twod().alpha*m2.get_twod().alpha;
      for (size_t r=0; r<nr; r++) {
        float *m3p=m3_d+r*m3_rs;
        for (size_t c=0; c<nc; c++) {
          const float *m1p=m1_d+r*m1_rs;
          const float *m2p=m2_d+c*m2_cs;
          float tot=0.0;
          for (size_t i=0; i<k; i++) {
            tot += *m1p * *m2p;
            m1p+=m1_cs;
            m2p+=m2_rs;
          }
          *m3p+=tot*fac;
          m3p+=m3_cs;
        }
      }
    } else {
      /* use cblas */
      const float ONE=float(1.0);
      const float alpha=(m1.get_twod().alpha * m2.get_twod().alpha);
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nr, nc, k, alpha, m1.get_twod().twod.data(), m1.nrows(), m2.get_twod().twod.data(), m2.nrows(), ONE, m3.data(), m3.nrows());
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

}

#else
////
// Implemententation of  matrix specializations without atlas/blais functionality
////

namespace mtl {
  
  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const double *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=1.0*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const double *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const double *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=m2.ncols();
    size_t m2_cs=1;
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*1.0;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    double *m3_d=m3.data();
    double fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      double *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        double tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=m3.ncols();
    size_t m3_cs=1;
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m2, column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m3)
  {
    size_t nr=m1.nrows();
    size_t nc=m2.ncols();
    size_t k=m1.ncols();
    assert(k==m2.nrows());
    assert(nr==m3.nrows());
    assert(nc==m3.ncols());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    size_t m2_rs=1;
    size_t m2_cs=m2.nrows();
    const float *m2_d=m2.get_twod().twod.data();
    size_t m3_rs=1;
    size_t m3_cs=m3.nrows();
    float *m3_d=m3.data();
    float fac=m1.get_twod().alpha*m2.get_twod().alpha;
    for (size_t r=0; r<nr; r++) {
      float *m3p=m3_d+r*m3_rs;
      for (size_t c=0; c<nc; c++) {
        const float *m1p=m1_d+r*m1_rs;
        const float *m2p=m2_d+c*m2_cs;
        float tot=0.0;
        for (size_t i=0; i<k; i++) {
          tot += *m1p * *m2p;
          m1p+=m1_cs;
          m2p+=m2_rs;
        }
        *m3p+=tot*fac;
        m3p+=m3_cs;
      }
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<double, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef double E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const double *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const double *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const row_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<row_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=m1.ncols();
    size_t m1_cs=1;
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=1.0*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> > &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=1.0*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< double >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const double *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const double *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*1.0;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float > &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*1.0;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< double >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    double *v3_d=v3.rep.data();
    double fac=m1.get_twod().alpha*v2.scale;
    double *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      double tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float > &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

  void mult(const column_matrix<gen_dense2D<float, gen_rect_offset<0, 0>, 0, 0 >, gen_rect_indexer<column_orien, 0, 0, size_t> >::scaled_type &m1, const dense1D< float >::scaled_type &v2, dense1D< float >::scaled_type &v3)
  {
    typedef float E;
    size_t nr=m1.nrows();
    size_t nc=m1.ncols();
    assert(nr==v3.size());
    assert(nc==v2.size());
    size_t m1_rs=1;
    size_t m1_cs=m1.nrows();
    const float *m1_d=m1.get_twod().twod.data();
    const float *v2_d=v2.rep.data();
    float *v3_d=v3.rep.data();
    float fac=m1.get_twod().alpha*v2.scale;
    float *v3p=v3_d;
    for (size_t r=0; r<nr; r++) {
      const float *m1p=m1_d+r*m1_rs;
      const float *v2p=v2_d;
      float tot=0.0;
      for (size_t c=0; c<nc; c++) {
        tot += *m1p * *v2p;
        m1p+=m1_cs;
        v2p++;
      }
      *v3p=tot*fac;
      v3p++;
    }
  }

}

#endif
