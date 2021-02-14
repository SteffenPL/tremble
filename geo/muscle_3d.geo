// Gmsh project created on Sun Apr  8 11:35:30 2018
//+
lc = DefineNumber[ 0.1, Name "Parameters/lc" ];
//+
R_1 = DefineNumber[ 0.3, Name "Parameters/R_1" ];
R_2 = DefineNumber[ 0.5, Name "Parameters/R_2" ];
R_3 = DefineNumber[ 0.2, Name "Parameters/R_3" ];
h = DefineNumber[ 1., Name "Parameters/h" ];
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {-R_1*0.3, R_1*0.3, 0, lc};
//+
Point(3) = {-R_1*0.3, R_1*0.9, 0, lc};
//+
Point(4) = {R_1*0.3, R_1*1.1, 0, lc};
//+
Point(5) = {R_1*0.7, R_1*0.8, 0, lc};
//+
Point(6) = {R_1*0.6, R_1*0.5, 0, lc};
//+
Point(7) = {R_1*0.3, R_1*0.4, 0, lc};
//+
Point(8) = {R_1*0.3, R_1*0.2, 0, lc};
//+
Point(10) = {-0.1, 0.1, 0.5*h, lc};
//+
Point(11) = {-0.35, 0.4, 0.5*h, lc};
//+
Point(12) = {-0.4, 0.75, 0.5*h, lc};
//+
Point(13) = {-0.25, 1, 0.5*h, lc};
//+
Point(14) = {0.15, 1.15, 0.5*h, lc};
//+
Point(15) = {0.55, 1.05, 0.5*h, lc};
//+
Point(16) = {0.75, 0.8, 0.5*h, lc};
//+
Point(17) = {0.6, 0.5, 0.5*h, lc};
//+
Point(18) = {0.35, 0.35, 0.5*h, lc};
//+
Point(19) = {0.2, 0.1, 0.5*h, lc};
//+
Point(20) = {0.05, 0.1, h, lc};
//+
Point(21) = {-0.2, 0.3, h, lc};
//+
Point(22) = {-0.2, 0.65, h, lc};
//+
Point(23) = {0.05, 0.9, h, lc};
//+
Point(24) = {0.35, 1.05, h, lc};
//+
Point(25) = {0.6, 1, h, lc};
//+
Point(26) = {0.75, 0.85, h, lc};
//+
Point(27) = {0.6, 0.6, h, lc};
//+
Point(28) = {0.35, 0.5, h, lc};
