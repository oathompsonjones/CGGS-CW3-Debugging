#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <queue>
#include "readOFF.h"



using namespace Eigen;
using namespace std;

bool isAnimating = false;
double timeStep = 0.02;
int numSpheres = 100;
double CRCoeff = 0.0;
double radius = 0.1;
double gravityConstant = 0.1;
MatrixXd spherePoses, sphereVelocities;
MatrixXd currV, origV;
MatrixXi F, fullF;


polyscope::SurfaceMesh* psMesh;




void callback_function() {
    ImGui::PushItemWidth(50);
    
    ImGui::TextUnformatted("Is Animating");
    ImGui::Separator();
    bool changed = ImGui::Checkbox("Is Animating", &isAnimating);
    ImGui::PopItemWidth();
    if (!isAnimating)
        return;
    
    //Collecting all gravitational forces
    //Earth gravity
    MatrixXd forces = MatrixXd::Zero(spherePoses.rows(),3);
    forces.col(1).array() = -9.8;
    //Pairwise gravitational forces
    double sqrRadius = radius*radius;
    for (int i=0;i<numSpheres;i++){
        for (int j=i+1;j<numSpheres;j++){
            RowVector3d forceDirection = spherePoses.row(j) - spherePoses.row(i);
            double sqrDistance = forceDirection.squaredNorm();
            if (sqrDistance > 4.1*sqrRadius){  //Otherwise the objects are too close and it will overshoot (tunnelling)
                forceDirection.normalize();
                forces.row(i).array() = forces.row(i) - forceDirection * gravityConstant / sqrDistance;
                forces.row(j).array() = forces.row(j) + forceDirection * gravityConstant / sqrDistance;
            }
        }
    }
    
    //Semi-implicit integration
    sphereVelocities += timeStep * forces;
    spherePoses += timeStep * sphereVelocities;
    
    //Doing collision detection and resolution
    for (int i=0;i<numSpheres;i++){
        for (int j=i+1;j<numSpheres;j++){
            RowVector3d colNormal = spherePoses.row(j) - spherePoses.row(i);
            double sqrDistance = colNormal.squaredNorm();
            //Explicitly finding out intersection between any two spheres
            if (sqrDistance < 2 * sqrRadius){
                colNormal.normalize();
                double colDist = 2.0 * radius - sqrt(sqrDistance);
                //Resolving interpenetration equally since masses and radii are equal
                spherePoses.row(i) = spherePoses.row(i) - colNormal * colDist / 2.0;
                spherePoses.row(j) = spherePoses.row(j) + colNormal * colDist / 2.0;
                //Hardcoding the linear velocity resolution (without inertia tensor) for mass = 1kg.
                double velBefore = (sphereVelocities.row(j) - sphereVelocities.row(i)).dot(colNormal);
                RowVector3d jn = -(1+CRCoeff) * velBefore*colNormal/2.0;
                sphereVelocities.row(i) -= jn;
                sphereVelocities.row(j) += jn;
            }
        }
        
        //Resolving collision between each sphere and the ground
        if (spherePoses(i,1) < radius){
            spherePoses(i, 1) = radius;
            if (sphereVelocities(i, 1) > 0.0)
                sphereVelocities(i, 1) = - CRCoeff * sphereVelocities(i, 1);
        }
    }
    
    for (int i=0;i<numSpheres;i++){
        currV.block(i*origV.rows(),0,origV.rows(), 3) = origV.rowwise() + spherePoses.row(i);
        fullF.block(i*F.rows(), 0, F.rows(),3) = F.array() + origV.rows()*i;
    }
    
    psMesh->updateVertexPositions(currV);
}


int main()
{
    readOFF(DATA_PATH "/spherers.off",origV, F);
    origV *= 2.0*radius;  //This is just for GUI; you don't need to touch it.
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    // Create a matrix for sphere positions (numSpheres x 3)
    spherePoses.resize(numSpheres, 3);

    // Fill the matrix with random values
    for (int i = 0; i < numSpheres; ++i) {
        for (int j = 0; j < 3; ++j) {
            spherePoses(i, j) = dist(gen);
        }
    }

    spherePoses.col(1).array() += 2.2;
    sphereVelocities = MatrixXd::Zero(spherePoses.rows(),3);
    
    //GUI stuff
    currV = MatrixXd::Zero(numSpheres*origV.rows(), 3);
    fullF.resize(numSpheres*F.rows(),3);
    for (int i=0;i<numSpheres;i++){
        currV.block(i*origV.rows(),0,origV.rows(), 3) = origV.rowwise() + spherePoses.row(i);
        fullF.block(i*F.rows(), 0, F.rows(),3) = F.array() + origV.rows()*i;
    }

    
    polyscope::init();
    psMesh = polyscope::registerSurfaceMesh("Balls", currV, fullF);
    psMesh->setSmoothShade(true);
    polyscope::state::lengthScale = 1.0;
    polyscope::state::boundingBox =
    std::tuple<glm::vec3, glm::vec3>{ {-5.0, 0, -5.0}, {5.0, 5.0, 5.0}};
    polyscope::options::groundPlaneHeightFactor = 0.0;
    
    
    
    polyscope::state::userCallback = callback_function;
    
    polyscope::show();
    
}

