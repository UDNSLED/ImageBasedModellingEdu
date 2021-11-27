//
// Created by caoqi on 2018/8/31.
//


//3D:  1.36939, -1.17123, 7.04869
//obs: 0.180123 -0.156584


#include "sfm/bundle_adjustment.h"
/*
 * This function computes the Jacobian entries for the given camera and
 * 3D point pair that leads to one observation.
 *
 * The camera block 'cam_x_ptr' and 'cam_y_ptr' is:
 * - ID 0: Derivative of focal length f
 * - ID 1-2: Derivative of distortion parameters k0, k1
 * - ID 3-5: Derivative of translation t0, t1, t2
 * - ID 6-8: Derivative of rotation w0, w1, w2
 *
 * The 3D point block 'point_x_ptr' and 'point_y_ptr' is:
 * - ID 0-2: Derivative in x, y, and z direction.
 *
 * The function that leads to the observation is given as follows:
 *
 *   u = f * D(x,y) * x  (image observation x coordinate)
 *   v = f * D(x,y) * y  (image observation y coordinate)
 *
 * with the following definitions:
 *
 *   xc = R0 * X + t0  (homogeneous projection)
 *   yc = R1 * X + t1  (homogeneous projection)
 *   zc = R2 * X + t2  (homogeneous projection)
 *   x = xc / zc  (central projection)
 *   y = yc / zc  (central projection)
 *   D(x, y) = 1 + k0 (x^2 + y^2) + k1 (x^2 + y^2)^2  (distortion)
 */

 /**
  * /description 给定一个相机参数和一个三维点坐标，求解雅各比矩阵，即公式中的df(theta)/dtheta
  * @param cam       相机参数
  * @param point     三维点坐标
  * @param cam_x_ptr 重投影坐标x 相对于相机参数的偏导数，相机有9个参数： [0] 焦距f; [1-2] 径向畸变系数k1, k2; [3-5] 平移向量 t1, t2, t3
  *                                                               [6-8] 旋转矩阵（角轴向量）
  * @param cam_y_ptr    重投影坐标y 相对于相机参数的偏导数，相机有9个参数
  * @param point_x_ptr  重投影坐标x 相对于三维点坐标的偏导数
  * @param point_y_ptr  重投影坐标y 相对于三维点坐标的偏导数
  */

void jacobian(sfm::ba::Camera const& cam,
              sfm::ba::Point3D const& point,
              double* cam_x_ptr, double* cam_y_ptr,
              double* point_x_ptr, double* point_y_ptr)
{
    const double f = cam.focal_length;
    const double *R = cam.rotation;
    const double *t = cam.translation;
    const double *X = point.pos;
    const double k0 = cam.distortion[0];
    const double k1 = cam.distortion[1];

    const double xc = R[0] * X[0] + R[1] * X[1] + R[2] * X[2] + t[0];
    const double yc = R[3] * X[0] + R[4] * X[1] + R[5] * X[2] + t[1];
    const double zc = R[6] * X[0] + R[7] * X[1] + R[8] * X[2] + t[2];

    const double x = xc / zc;
    const double y = yc / zc;
    const double r2 = x * x + y * y;
    const double d = 1 + (k0 + k1 * r2) * r2;


    // 相机焦距的偏导数
    cam_x_ptr[0] = d * x;
    cam_y_ptr[0] = d * y;

    // 相机径向畸变的偏导数
    cam_x_ptr[1] = f * x * r2;
    cam_x_ptr[2] = f * x * r2 * r2;
    cam_y_ptr[1] = f * y * r2;
    cam_y_ptr[2] = f * y * r2 * r2;

    // 相机平移向量的偏导数
    double du_dd = f * x;
    double dv_dd = f * y;

    double dd_dxc =  (k0 + 2 * k1 * r2) * 2 * x / zc;
    double dd_dyc =  (k0 + 2 * k1 * r2) * 2 * y / zc;
    double dd_dzc = -(k0 + 2 * k1 * r2) * 2 * r2 / zc;
    
    double du_dx = f * d;
    double dv_dy = f * d;

    double dx_dxc = 1 / zc;
    double dx_dyc = 0;
    double dx_dzc = -x / zc;
    double dy_dxc = 0;
    double dy_dyc = 1/zc;
    double dy_dzc = -y / zc;

    double du_dxc = du_dd * dd_dxc + du_dx * dx_dxc;
    double du_dyc = du_dd * dd_dyc + du_dx * dx_dyc;
    double du_dzc = du_dd * dd_dzc + du_dx * dx_dzc;

    double dv_dxc = dv_dd * dd_dxc + dv_dy * dy_dxc;
    double dv_dyc = dv_dd * dd_dyc + dv_dy * dy_dyc;
    double dv_dzc = dv_dd * dd_dzc + dv_dy * dy_dzc;

    cam_x_ptr[3] = du_dxc;
    cam_x_ptr[4] = du_dyc;
    cam_x_ptr[5] = du_dzc;
    cam_y_ptr[3] = dv_dxc;
    cam_y_ptr[4] = dv_dyc;
    cam_y_ptr[5] = dv_dzc;

    double dxc_dw0 = 0;
    double dxc_dw1 = (R[6] * X[0] + R[7] * X[1] + R[8] * X[2]);
    double dxc_dw2 = -(R[3] * X[0] + R[4] * X[1] + R[5] * X[2]);

    double dyc_dw0 = - (R[6] * X[0] + R[7] * X[1] + R[8] * X[2]);
    double dyc_dw1 = 0;
    double dyc_dw2 = R[0] * X[0] + R[1] * X[1] + R[2] * X[2];

    double dzc_dw0 = (R[3] * X[0] + R[4] * X[1] + R[5] * X[2]);
    double dzc_dw1 = -(R[0] * X[0] + R[1] * X[1] + R[2] * X[2]);
    double dzc_dw2 = 0;

    double du_dw0 = du_dyc * dyc_dw0 + du_dzc * dzc_dw0;
    double du_dw1 = du_dxc * dxc_dw1 + du_dzc * dzc_dw1;
    double du_dw2 = du_dxc * dxc_dw2 + du_dyc * dyc_dw2;

    double dv_dw0 = dv_dyc * dyc_dw0 + dv_dzc * dzc_dw0;
    double dv_dw1 = dv_dxc * dxc_dw1 + dv_dzc * dzc_dw1;
    double dv_dw2 = dv_dxc * dxc_dw2 + dv_dyc * dyc_dw2;

    // 相机旋转矩阵的偏导数
    cam_x_ptr[6] = du_dw0;
    cam_x_ptr[7] = du_dw1;
    cam_x_ptr[8] = du_dw2;
    cam_y_ptr[6] = dv_dw0;
    cam_y_ptr[7] = dv_dw1;
    cam_y_ptr[8] = dv_dw2;

    double dxc_dX = R[0];
    double dxc_dY = R[1];
    double dxc_dZ = R[2];

    double dyc_dX = R[3];
    double dyc_dY = R[4];
    double dyc_dZ = R[5];

    double dzc_dX = R[6];
    double dzc_dY = R[7];
    double dzc_dZ = R[8];

    // 三维点的偏导数
    point_x_ptr[0] = du_dxc * dxc_dX + du_dyc * dyc_dX + du_dzc * dzc_dX;
    point_x_ptr[1] = du_dxc * dxc_dY + du_dyc * dyc_dY + du_dzc * dzc_dY;
    point_x_ptr[2] = du_dxc * dxc_dZ + du_dyc * dyc_dZ + du_dzc * dzc_dZ;
    point_y_ptr[0] = dv_dxc * dxc_dX + dv_dyc * dyc_dX + dv_dzc * dzc_dX;
    point_y_ptr[1] = dv_dxc * dxc_dY + dv_dyc * dyc_dY + dv_dzc * dzc_dY;
    point_y_ptr[2] = dv_dxc * dxc_dZ + dv_dyc * dyc_dZ + dv_dzc * dzc_dZ;

}
int main(int argc, char*argv[])
{

    sfm::ba::Camera cam;
    cam.focal_length  =  0.919654;
    cam.distortion[0] = -0.108298;
    cam.distortion[1] =  0.103775;

    cam.rotation[0] = 0.999999;
    cam.rotation[1] = -0.000676196;
    cam.rotation[2] = -0.0013484;
    cam.rotation[3] = 0.000663243;
    cam.rotation[4] = 0.999949;
    cam.rotation[5] = -0.0104095;
    cam.rotation[6] = 0.00135482;
    cam.rotation[7] = 0.0104087;
    cam.rotation[8] = 0.999949;

    cam.translation[0]=0.00278292;
    cam.translation[1]=0.0587996;
    cam.translation[2]=-0.127624;

    sfm::ba::Point3D pt3D;
    pt3D.pos[0]= 1.36939;
    pt3D.pos[1]= -1.17123;
    pt3D.pos[2]= 7.04869;

    double cam_x_ptr[9]={0};
    double cam_y_ptr[9]={0};
    double point_x_ptr[3]={0};
    double point_y_ptr[3]={0};

    jacobian(cam, pt3D, cam_x_ptr, cam_y_ptr, point_x_ptr, point_y_ptr);


   std::cout<<"Result is :"<<std::endl;
    std::cout<<"cam_x_ptr: ";
    for(int i=0; i<9; i++){
        std::cout<<cam_x_ptr[i]<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"cam_y_ptr: ";
    for(int i=0; i<9; i++){

        std::cout<<cam_y_ptr[i]<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"point_x_ptr: ";
    std::cout<<point_x_ptr[0]<<" "<<point_x_ptr[1]<<" "<<point_x_ptr[2]<<std::endl;

    std::cout<<"point_y_ptr: ";
    std::cout<<point_y_ptr[0]<<" "<<point_y_ptr[1]<<" "<<point_y_ptr[2]<<std::endl;


    std::cout<<"\nResult should be :\n"
       <<"cam_x_ptr: 0.195942 0.0123983 0.000847141 0.131188 0.000847456 -0.0257388 0.0260453 0.95832 0.164303\n"
       <<"cam_y_ptr: -0.170272 -0.010774 -0.000736159 0.000847456 0.131426 0.0223669 -0.952795 -0.0244697 0.179883\n"
       <<"point_x_ptr: 0.131153 0.000490796 -0.0259232\n"
       <<"point_y_ptr: 0.000964926 0.131652 0.0209965\n";


    return 0;
}
