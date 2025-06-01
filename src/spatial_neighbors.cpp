#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// helper to generate offsets based on connectivity
static std::vector< std::array<int,3> >
make_offsets(int connectivity) {
  std::vector< std::array<int,3> > offs;
  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        int manhattan = std::abs(dx) + std::abs(dy) + std::abs(dz);
        if (connectivity == 6) {
          if (manhattan == 1) offs.push_back({{dx,dy,dz}});
        } else if (connectivity == 18) {
          if (manhattan <= 2) offs.push_back({{dx,dy,dz}});
        } else if (connectivity == 26) {
          offs.push_back({{dx,dy,dz}});
        } else {
          stop("Unsupported connectivity");
        }
      }
    }
  }
  return offs;
}

// [[Rcpp::export]]
List spatial_neighbors_cpp(const LogicalVector& mask,
                           const IntegerVector& dims,
                           int connectivity) {
  if (dims.size() != 3)
    stop("dims must have length 3");
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
  int total = nx * ny * nz;
  if (mask.size() != total)
    stop("mask length does not match dims");

  // mapping from linear index to voxel id
  std::vector<int> mapping(total, 0);
  std::vector< std::array<int,3> > coords_vec;
  coords_vec.reserve(total);

  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        int lin = x + y*nx + z*nx*ny;
        if (mask[lin]) {
          mapping[lin] = coords_vec.size() + 1; // 1-based
          coords_vec.push_back({{x+1,y+1,z+1}});
        }
      }
    }
  }

  int n_vox = coords_vec.size();
  List neighbors(n_vox);
  std::vector<int> i_vec;
  std::vector<int> j_vec;
  std::vector<double> x_vec;
  i_vec.reserve(n_vox * 7);
  j_vec.reserve(n_vox * 7);
  x_vec.reserve(n_vox * 7);

  std::vector< std::array<int,3> > offs = make_offsets(connectivity);

  for (int v = 0; v < n_vox; ++v) {
    int x = coords_vec[v][0];
    int y = coords_vec[v][1];
    int z = coords_vec[v][2];
    std::vector<int> nb_local;
    nb_local.reserve(offs.size());
    for (auto& o : offs) {
      int nxp = x + o[0];
      int nyp = y + o[1];
      int nzp = z + o[2];
      if (nxp >= 1 && nxp <= nx &&
          nyp >= 1 && nyp <= ny &&
          nzp >= 1 && nzp <= nz) {
        int lin = (nxp-1) + (nyp-1)*nx + (nzp-1)*nx*ny;
        int nb_idx = mapping[lin];
        if (nb_idx > 0) {
          nb_local.push_back(nb_idx);
        }
      }
    }
    IntegerVector nb(nb_local.begin(), nb_local.end());
    neighbors[v] = nb;

    int degree = nb_local.size();
    i_vec.push_back(v + 1);
    j_vec.push_back(v + 1);
    x_vec.push_back(static_cast<double>(degree));
    for (int nb_idx : nb_local) {
      i_vec.push_back(v + 1);
      j_vec.push_back(nb_idx);
      x_vec.push_back(-1.0);
    }
  }

  IntegerMatrix coords_out(n_vox, 3);
  for (int v = 0; v < n_vox; ++v) {
    coords_out(v,0) = coords_vec[v][0];
    coords_out(v,1) = coords_vec[v][1];
    coords_out(v,2) = coords_vec[v][2];
  }

  IntegerVector I(i_vec.begin(), i_vec.end());
  IntegerVector J(j_vec.begin(), j_vec.end());
  NumericVector X(x_vec.begin(), x_vec.end());

  return List::create(Named("neighbors") = neighbors,
                      Named("coords") = coords_out,
                      Named("i") = I,
                      Named("j") = J,
                      Named("x") = X);
}

// [[Rcpp::export]]
List laplacian_from_neighbors_cpp(List neighbors) {
  int n_voxels = neighbors.size();
  std::vector<int> i_vec;
  std::vector<int> j_vec;
  std::vector<double> x_vec;
  for (int v = 0; v < n_voxels; ++v) {
    IntegerVector nb = neighbors[v];
    int degree = nb.size();
    i_vec.push_back(v + 1);
    j_vec.push_back(v + 1);
    x_vec.push_back(static_cast<double>(degree));
    for (int k = 0; k < degree; ++k) {
      i_vec.push_back(v + 1);
      j_vec.push_back(nb[k]);
      x_vec.push_back(-1.0);
    }
  }
  IntegerVector I(i_vec.begin(), i_vec.end());
  IntegerVector J(j_vec.begin(), j_vec.end());
  NumericVector X(x_vec.begin(), x_vec.end());
  return List::create(Named("i")=I, Named("j")=J, Named("x")=X);
}

