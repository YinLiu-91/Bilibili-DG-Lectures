% main.m
% 参考文章：
%1. https://zhuanlan.zhihu.com/p/605184617
%2. https://zhuanlan.zhihu.com/p/604524727
clear;clc

global Nx dimPk NumGLP
Nx = 100;           % 单元数量100个
k = 2;              % 多项式空间Pk的次数，那么每个单元上需要求解3x3的矩阵
NumGLP = 5;         % Gauss-Legendre点的个数
dimPk = k + 1;
CFL = 0.2;          % CFL数

get_GLP             %  获取G-L点和权重

init_data           % 赋初值

get_basis           % 储存多项式基的函数值、导数值、左右端点值，以及质量

L2_Pro              % 用L2投影获取初值

RK3                 % 用3阶RK进行迭代

calculate_L2_Error  %计算数值解的L2误差
