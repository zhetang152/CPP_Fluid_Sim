#include"vector3D.h"
#include"grid.h"
#include"MACgrid.h"

void setupScene(MACGrid& grid) {
    auto& types = grid.celltypes();
    int nx = types.getWidth();
    int ny = types.getHeight();
    int nz = types.getDepth();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (i == 0 || i == nx - 1 || 
                    j == 0 || j == ny - 1 || 
                    k == 0 || k == nz - 1) {
                    types(i, j, k) = CellType::SOLID;
                } else {
                    types(i, j, k) = CellType::FLUID;
                }
            }
        }
    }
}
Grid<float> discrete_divergence(const MACGrid& grid){
    int nx = grid.celltypes().getWidth();
    int ny = grid.celltypes().getHeight();
    int nz = grid.celltypes().getDepth();
    float dx = grid.getDx();
    Grid<float> negetivedivergence(nx,ny,nz,0.0f);
    const auto& u = grid.u();
    const auto& v = grid.v();
    const auto& w = grid.w();
    const auto& celltypes = grid.celltypes();
    for(int k = 0; k<nz;++k){
        for(int j = 0; j<ny;++j){
            for(int i = 0;i<nx;++i){
                if(celltypes(i,j,k)==CellType::FLUID){
                    float u_right = u(i+1,j,k);
                    float u_left = u(i,j,k);
                    float v_top = v(i,j+1,k);
                    float v_bottom = v(i,j,k);
                    float w_front = w(i,j,k+1);
                    float w_back = w(i,j,k);
                    float divergence = (u_right-u_left)+(v_top-v_bottom)+(w_front-w_back);
                    negetivedivergence(i,j,k)=-divergence/dx;
                }
            }
        }
    }
    return negetivedivergence;
}