// Gmsh project created on Thu May 17 11:30:56 2018
//+
lc = DefineNumber[ 0.1, Name "Parameters/lc" ];
//+
top = DefineNumber[ 1, Name "Parameters/top" ];
//+
bot = DefineNumber[ 0, Name "Parameters/bot" ];

r1 = DefineNumber[ 0.3, Name "Parameters/r1" ];
r2 = DefineNumber[ 0.5, Name "Parameters/r2" ];
r3 = DefineNumber[ 0.3, Name "Parameters/r3" ];
//+
Point(1) = {r1*0, -r1*0.3, bot, lc};
//+
Point(2) = {r1*0.3, -r1*0.2, bot, lc};
//+
Point(3) = {r1*0.4, r1*0.1, bot, lc};
//+
Point(4) = {r1*0, r1*0.3, bot, lc};
//+
Point(5) = {-r1*0.4, r1*0.2, bot, lc};
//+
Point(6) = {-r1*0.3, -r1*0.2, bot, lc};


Point(7) = {r2*0, -r2*0.3, 0.5*bot+0.5*top, lc};
//+
Point(8) = {r2*0.3, -r2*0.2, 0.5*bot+0.5*top, lc};
//+
Point(9) = {r2*0.3, r2*0.1, 0.5*bot+0.5*top, lc};
//+
Point(10) = {r2*0, r2*0.4, 0.5*bot+0.5*top, lc};
//+
Point(11) = {-r2*0.3, r2*0.2, 0.5*bot+0.5*top, lc};
//+
Point(12) = {-r2*0.3, -r2*0.2, 0.5*bot+0.5*top, lc};


Point(13) = {r3*0, -r3*0.3, top, lc};
//+
Point(14) = {r3*0.3, -r3*0.2, top, lc};
//+
Point(15) = {r3*0.4, r3*0.1, top, lc};
//+
Point(16) = {r3*0, r3*0.3, top, lc};
//+
Point(17) = {-r3*0.4, r3*0.2, top, lc};
//+
Point(18) = {-r3*0.3, -r3*0.2, top, lc};
//+
Spline(1) = {3, 4, 5, 6, 1, 2, 3};
//+
Spline(2) = {9, 10, 11, 12, 7, 8, 9};
//+
Spline(3) = {15, 16, 17, 18, 13, 14, 15};
//+
Spline(4) = {4, 10, 16};
//+
Spline(5) = {3, 9, 15};
//+
Spline(6) = {5, 11, 17};
//+
Spline(7) = {18, 12, 6};
//+
Spline(8) = {1, 7, 13};
//+
Spline(9) = {14, 8, 2};
//+
SetFactory("Built-in");
//+
Coherence;
//+
SetFactory("OpenCASCADE");
Wire(1) = {8};
Extrude { Curve{1}; } Using Wire {1}

//+
SetFactory("OpenCASCADE");
//+
Wire(1) = {7};
Extrude { Curve{1}; } Using Wire {1}

//+
Wire(1) = {8};
Extrude { Curve{1}; Curve{2}; Curve{3}; } Using Wire {1}

//+
Extrude {0, 0, 1} {
  Curve{1}; Curve{8}; 
}
//+
Extrude {0, 0, 1} {
  Curve{1}; 
}
//+
Extrude {0, 0, 1} {
  Curve{8}; 
}
