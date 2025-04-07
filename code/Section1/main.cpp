#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <queue>
#include "readOFF.h"
#include "create_edge_list.h"
#include "compute_laplacian.h"
#include "slice_columns_sparse.h"
#include "set_diff.h"


using namespace Eigen;
using namespace std;

MatrixXi F, E, EF;
VectorXi boundEMask, boundVMask, boundVertices;
MatrixXd V, Hn;
SparseMatrix<double> d0, W;
VectorXd vorAreas, H;
vector<vector<int>> oneRings;  //every vector<int> is a single one ring (not including the original vertex)

MatrixXd ARAP_deformation(const MatrixXd& origV,
                          const MatrixXi& E,
                          const SparseMatrix<double>& d0,
                          const SparseMatrix<double>& W,
                          const vector<vector<int>>& oneRings,
                          const VectorXi& constHandles,
                          const MatrixXd& constPositions,
                          const int numIterations){
    
    vector<Matrix3d> R(origV.rows());
    for (int i=0;i<origV.rows();i++)
        R[i] = Matrix3d::Identity();
    
    MatrixXd g(E.rows(),3);
    
    //Constructing the fixed left-hand side.
    SparseMatrix<double> d0B = slice_columns_sparse(d0, constHandles);
    VectorXi allVertices(d0.cols());
    for (int i=0;i<d0.cols();i++)
        allVertices(i)=i;
    VectorXi varVertices = set_diff(allVertices, constHandles);
    SparseMatrix<double> d0I = slice_columns_sparse(d0, varVertices);
    SparseMatrix<double> A = d0I.transpose()*W*d0I;
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A);
    MatrixXd currV(origV.rows(), 3);
    for (int i=0;i<numIterations;i++){
        //Global step: rotate each edge ij by (Ri+Rj)/2 for the previous or initial iteration's R (see Lecture notes) and try to integrate to it
        for (int j=0;j<E.rows();j++){
            Matrix3d REdge = (R[E(j,0)]+R[E(j,1)])/2;
            g.row(j) = (origV.row(E(j,0)) - origV.row(E(j,1)))*REdge;
        }
        
        MatrixXd rhs = d0I.transpose() * (g - d0B*constPositions);
        MatrixXd x = solver.solve(rhs);
        for (int j=0;j<constPositions.rows();j++)
            currV.row(constHandles(j))=constPositions.row(j);
        for (int j=0;j<varVertices.size();j++)
            currV.row(varVertices(j))=x.row(j);
        
        //Local step: for existing currVertices and original positions origVertices, find the best fitting local one-ring-based rotation matrices R
        for (int j=0;j<origV.rows();j++){
            MatrixXd P(oneRings[j].size(),3), Q(oneRings[j].size(),3);
            Matrix3d S;
            for (int k=0;k<oneRings[j].size();k++){
                P.row(k) = origV.row(oneRings[j][k]) - origV.row(j);
                P.row(k) = origV.row(oneRings[j][k]) - origV.row(j);
                Q.row(k) = currV.row(oneRings[j][k]) - currV.row(j);
                Q.row(k) = currV.row(oneRings[j][k]) - currV.row(j);
                S = Q.transpose()*P;
                Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Matrix3d U = svd.matrixU();
                Eigen::Vector3d Sigma = svd.singularValues();
                Eigen::Matrix3d Vt = svd.matrixV().transpose();
                
                Matrix3d currR = U*Vt;
                if (currR.determinant()<0.0){
                    //check where the smallest singular values falls
                    int minValue, minIndex;
                    minValue = Sigma.minCoeff(&minIndex);
                    Matrix3d newSigma = Matrix3d::Identity();
                    newSigma(minIndex, minIndex) = -1;
                    currR = U*newSigma*Vt;
                }
                R[j] = currR;
            }
        }
    }
    return currV;
}
    


int main()
{
    readOFF(DATA_PATH "/horsers.off",V, F);
    create_edge_list(F, E, EF, boundEMask, boundVMask, boundVertices, oneRings);
    compute_laplacian(V, F, E, EF, boundEMask, d0, W, vorAreas);

    polyscope::init();
    polyscope::SurfaceMesh* psInputMesh = polyscope::registerSurfaceMesh("Input Mesh", V, F);
    polyscope::SurfaceMesh* psOutputMesh = polyscope::registerSurfaceMesh("Output Mesh", V, F);
    
    
    VectorXi constHandles(5);
    constHandles<<535, 780, 793, 827, 769;
    int numIterations = 25;

    MatrixXd constPositions(5,3);
    constPositions<<V.row(constHandles(0))+RowVector3d({0.01, 0.1, 0.01}),
                    V.row(constHandles(1))+RowVector3d({0, 0.09, -0.1}),
                    V.row(constHandles(2))+RowVector3d({0, 0.09, -0.11}),
                    V.row(constHandles(3)), V.row(constHandles(4));
    
    polyscope::registerPointCloud("const handles", constPositions);
    
    MatrixXd currVertices = ARAP_deformation(V, E, d0, W, oneRings, constHandles, constPositions, numIterations);

    psOutputMesh->updateVertexPositions(currVertices);
    
    polyscope::show();
    
}

