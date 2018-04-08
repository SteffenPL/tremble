// Gmsh project created on Sun Apr  8 12:10:40 2018
//+
T_1 = DefineNumber[ 1, Name "Parameters/T_1" ];
W_1 = DefineNumber[ 1, Name "Parameters/W_1" ];
H_1 = DefineNumber[ 0, Name "Parameters/H_1" ];

T_2 = DefineNumber[ 1, Name "Parameters/T_2" ];
W_2 = DefineNumber[ 1, Name "Parameters/W_2" ];
H_2 = DefineNumber[ 0.5, Name "Parameters/H_2" ];

T_3 = DefineNumber[ 1, Name "Parameters/T_3" ];
W_3 = DefineNumber[ 1, Name "Parameters/W_3" ];
H_3 = DefineNumber[ 1., Name "Parameters/H_3" ];

lc = DefineNumber[ 0.5, Name "Parameters/lc" ];
f_1 = DefineNumber[ 1., Name "Parameters/f_1" ];
f_2 = DefineNumber[ 0.5, Name "Parameters/f_2" ];
f_3 = DefineNumber[ 1., Name "Parameters/f_3" ];

lc_1=lc*f_1;
lc_2=lc*f_2;
lc_3=lc*f_3;

For i In {1:3}
    Point(1+4*(i-1)) = {-0.5*W~{i},     0, H~{i}, lc~{i}};
    Point(2+4*(i-1)) = { 0.5*W~{i},     0, H~{i}, lc~{i}};
    Point(3+4*(i-1)) = { 0.5*W~{i}, T~{i}, H~{i}, lc~{i}};
    Point(4+4*(i-1)) = {-0.5*W~{i}, T~{i}, H~{i}, lc~{i}};
EndFor
//+
Line(1) = {1, 2};
//+
Line(2) = {5, 6};
//+
Line(3) = {9, 10};
//+
BSpline(4) = {1, 4, 3, 2};
//+
BSpline(5) = {5, 8, 7, 6};
//+
BSpline(6) = {9, 12, 11, 10};
//+
Line(7) = {1, 5};
//+
Line(8) = {2, 6};
//+
Line(9) = {6, 10};
//+
Line(10) = {5, 9};
//+
Line Loop(1) = {1, -4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {3, -6};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {10, 3, -9, -2};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {7, 2, -8, -1};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {4, 8, -5, -7};
//+
Surface(5) = {5};
//+
Line Loop(6) = {5, 9, -6, -10};
//+
Surface(6) = {6};
//+
Surface Loop(1) = {2, 3, 6, 5, 1, 4};
//+
Volume(1) = {1};
