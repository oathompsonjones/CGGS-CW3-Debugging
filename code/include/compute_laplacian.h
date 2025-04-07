#ifndef COMPUTE_LAPLACIAN_HEADER_FILE
#define COMPUTE_LAPLACIAN_HEADER_FILE

#include <Eigen/Dense>

void compute_laplacian(const Eigen::MatrixXd& V,
                       const Eigen::MatrixXi& F,
                       const Eigen::MatrixXi& E,
                       const Eigen::MatrixXi& EF,
                       const Eigen::VectorXi& boundEMask,
                       Eigen::SparseMatrix<double>& d0,
                       Eigen::SparseMatrix<double>& W,
                       Eigen::VectorXd& vorAreas){
    
    using namespace Eigen;
    using namespace std;
    d0.resize(E.rows(), V.rows());
    W.resize(E.rows(), E.rows());
    vector<Triplet<double>> d0Tris, WTris;
    for (int i=0;i<E.rows();i++){
        d0Tris.push_back(Triplet<double>(i, E(i,0),-1.0));
        d0Tris.push_back(Triplet<double>(i, E(i,1), 1.0));
        
        RowVector3d vk = V.row(F(EF(i,0), (EF(i,1)+1)%3));
        RowVector3d vi = V.row(F(EF(i,0), (EF(i,1)+2)%3));
        RowVector3d vj = V.row(F(EF(i,0), EF(i,1)));
        double cosj = (vk - vj).dot(vi - vj);
        double sinj = ((vk - vj).cross(vi - vj)).norm();
        WTris.push_back(Triplet<double>(i,i,0.5 * cosj / sinj));
        if (!boundEMask(i)){
            RowVector3d vl = V.row(F(EF(i,2), EF(i,3)));
            double cosl = (vi - vl).dot(vk - vl);
            double sinl = ((vi - vl).cross(vk - vl)).norm();
            WTris.push_back(Triplet<double>(i,i,0.5*cosl/sinl));
        }
    }
    d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
    W.setFromTriplets(WTris.begin(), WTris.end());
    
    //computing Voronoi areas
    vorAreas = VectorXd::Zero(V.rows());
    for (int i=0;i<F.rows();i++){
        RowVector3d v1 = V.row(F(i,0));
        RowVector3d v2 = V.row(F(i,1));
        RowVector3d v3 = V.row(F(i,2));
        double faceArea = ((v2-v1).cross(v3-v1)).norm()*0.5;
        
        for (int j=0;j<3;j++)
            vorAreas(F(i,j))+=faceArea / 3.0;
    }
}



#endif
