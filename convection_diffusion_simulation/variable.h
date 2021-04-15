#pragma once

const double lx = 1;
const double ly = 1;
const int nx = 12;
const int ny = 12;
const int n = nx * ny;
const double dx = lx / (nx - 2);
const double dy = ly / (ny - 2);

const double gamma = 0.1;   //扩散系数
const double rho = 1;       //密度
const double u = 1;       //速度

double D = gamma * (dx * 1) / dx;    //第一个dx面积，第二个距离
double F = rho * u * (dx * 1);       //这里dx代表面积

double dt = 0.01;
double ap0 = rho * 1 * dx * dy / dt;

