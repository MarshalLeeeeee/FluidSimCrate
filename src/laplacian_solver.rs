/// Module for sparse laplacian solver
use ndarray as nd;
use sprs_ldl::{Ldl, LdlNumeric};
use sprs::{CsMat, TriMat};

/// Laplacian solver that holds sparse laplacian operator as A, where the linear system is Ax=b
pub struct LaplacianSolver {
    system: LdlNumeric<f64, usize>,
}
impl LaplacianSolver {
    /// Create instance of laplacian solver with given width and height
    pub fn new(w: usize, h: usize) -> Self {
        let sz = w * h;
        let mut pa = TriMat::<f64>::new((sz, sz));
        for i in 0..w {
            for j in 0..h {
                let mut cnt = 0;
                if i > 0 {
                    pa.add_triplet(i*h+j, (i-1)*h+j, 1_f64);
                    cnt += 1;
                }
                if i < w-1 {
                    pa.add_triplet(i*h+j, (i+1)*h+j, 1_f64);
                    cnt += 1;
                }
                if j > 0 {
                    pa.add_triplet(i*h+j, i*h+j-1, 1_f64);
                    cnt += 1;
                }
                if j < h-1 {
                    pa.add_triplet(i*h+j, i*h+j+1, 1_f64);
                    cnt += 1;
                }
                pa.add_triplet(i*h+j, i*h+j, -cnt as f64);
            }
        }
        let pa: CsMat<_> = pa.to_csr();
        let ldl = Ldl::default();
        Self {
            system: ldl.numeric(pa.view()).unwrap()
        }
    }

    /// Solve Ax=b given the dense vector b
    pub fn solve(&self, b: nd::Array1::<f64>) -> nd::Array1::<f64> {
        self.system.solve(b)
    }
}