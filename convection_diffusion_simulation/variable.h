#pragma once

const double lx = 1;
const double ly = 1;
const int nx = 12;
const int ny = 12;
const int n = nx * ny;
const double dx = lx / (nx - 2);
const double dy = ly / (ny - 2);

const double gamma = 0.1;   //��ɢϵ��
const double rho = 1;       //�ܶ�
const double u = 1;       //�ٶ�

double D = gamma * (dx * 1) / dx;    //��һ��dx������ڶ�������
double F = rho * u * (dx * 1);       //����dx�������

double dt = 0.01;
double ap0 = rho * 1 * dx * dy / dt;

