#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;

/* Rectangle Format:
___________________   -y1
|        |        |
|    1   |    2   |
|--------|--------|   -y2
|    0   |    3   |
x0_______x2_______x1  -y0
*/
int const MINLEVEL = 3;
int const MAXLEVEL = 15;
int i = 0;
int j = 0; //counters for pow function
ostream *rSizePtr;

//The 0-value is the leftmost endpoint, the 1-value is the rightmost, and the 2-value is midpoint
typedef struct rectangle {
    double x[3];
    double y[3];
} rect;

//Returns area of rectangle r
double area(rect r) {
    return fabs(r.x[0] - r.x[1]) * fabs(r.y[0] - r.y[1]);
}

double power3D(double x, double y) {
    return pow(x, i) * pow(y, j);
}

double unitCircle(double x, double y) {
    return (x * x + y * y < 1) ? 1.0 : 0.0;
}

//Used for printing out good rectangle sizes
void outRect(rect r) {
    *rSizePtr << '\n' << r.x[0] << ' ' << r.y[0] << '\n' << r.x[0] << ' ' << r.y[1] << '\n' << r.x[1]
              << ' ' << r.y[1] << '\n' << r.x[1] << ' ' << r.y[0] << '\n' << r.x[0] << ' ' << r.y[0] << endl;
}

//Calculates fine approximation of integral
double fineInt(double (*f)(double, double), rect r, double *fvals) {
    double int1, int2, int3, int4;
    double fmid = f(r.x[2], r.y[2]);
    int1 = (fvals[0] + f(r.x[0], r.y[2]) + f(r.x[2], r.y[0]) + fmid);
    int2 = (fvals[1] + f(r.x[2],r.y[1]) + fmid + f(r.x[0], r.y[2]));
    int3 = (fvals[2] + f(r.x[2],r.y[1]) + fmid + f(r.x[1], r.y[2]));
    int4 = (fvals[3] + f(r.x[2], r.y[0]) + fmid + f(r.x[1], r.y[2]));
    return 0.25 * (int1 + int2 + int3 + int4) * (area(r) * 0.25);
}

//Does recursive calling of new rectangles and also splits up original into 4.
double tensorTrap(double (*f)(double, double), rect r, double tol, int lev) {
    double fvals[4];
    fvals[0] = f(r.x[0], r.y[0]);
    fvals[1] = f(r.x[0], r.y[1]);
    fvals[2] = f(r.x[1], r.y[1]);
    fvals[3] = f(r.x[1], r.y[0]);
    double coarseApprox = 0.25 * (fvals[0] + fvals[1] + fvals[2] + fvals[3]) * area(r);
    double fineApprox = fineInt(f, r, fvals);
    if(((fabs(fineApprox - coarseApprox) < 3 * tol) && lev >= MINLEVEL) || lev >= MAXLEVEL) {
        outRect(r);
        return fineApprox;
    }

    rect r0 = {{r.x[0], r.x[2], (r.x[0] + r.x[2]) / 2}, {r.y[0], r.y[2], (r.y[2] + r.y[0]) / 2}};
    rect r1 = {{r.x[0], r.x[2], (r.x[0] + r.x[2]) / 2}, {r.y[2], r.y[1], (r.y[2] + r.y[1]) / 2}};
    rect r2 = {{r.x[2], r.x[1], (r.x[1] + r.x[2]) / 2}, {r.y[2], r.y[1], (r.y[2] + r.y[1]) / 2}};
    rect r3 = {{r.x[2], r.x[1], (r.x[1] + r.x[2]) / 2,}, {r.y[0], r.y[2], (r.y[2] + r.y[0]) / 2}};

    return tensorTrap(f, r0, tol / 4, lev + 1) + tensorTrap(f, r1, tol / 4, lev + 1)
           + tensorTrap(f, r2, tol / 4, lev + 1) + tensorTrap(f, r3, tol / 4, lev + 1);
}

int main() {
    rect r = {{0,1,0.5}, {0,1,0.5}};
    for(i = 0; i <= 4; i ++) {
        for(j = i; j <= 4; j ++) {
            ofstream rectSizeOut(("x" + to_string(i) + "y" + to_string(j)).c_str());
            rSizePtr = &rectSizeOut;
            printf("x^%dy^%d: %.7f\n", i, j, tensorTrap(&power3D, r, 1/1000.0, 1));
            rectSizeOut.close();
        }
    }
    ofstream circleRect("circle");
    rSizePtr = &circleRect;
    rect r2 = {{-1, 1, 0}, {-1, 1, 0}};
    printf("Unit Circle: %.7f\n", tensorTrap(&unitCircle, r2, 1/1000.0, 1));
    system("gnuplot ./plot.gnu");
}