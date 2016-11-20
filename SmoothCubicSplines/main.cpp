#import "hermiteCubic.h"
#include <fstream>

using namespace std;

double maxError(hermiteCubic h) { //Approxs function using cspline to find max error
    double err = 0;
    for(int i = 0; i < 1000; i ++) {
        double x = i * M_PI_2 / 1000;
        err = max(fabs(h.evaluate(x) - sin(x)), err);
    }
    return err;
}

void test1(hermiteCubic h) {
    vector<double> answers = smoothCubic_DEC(h);
    ofstream t1("Part1Curve");
    for(int i = 0; i < h.n; i ++) {
        cout << "g'(" << h.p[i] << ") = " << answers[i] << endl;
    }
    for(int i = 0; i <= 25; i ++) {
        double s = i * 1.0 / 5;
        t1 << s << '\t' << h.evaluate(s) << '\n';
    }
}

void test2(hermiteCubic h) {
    vector<double> answers = not_a_knot(h);
    ofstream t2("Part2Curve");
    for(int i = 0; i < h.n; i ++) {
        cout << "g'(" << h.p[i] << ") = " << answers[i] << endl;
    }
    for(int i = 0; i <= 50; i ++) {
        double s = i * 1.0 / 10.4;
        t2 << s << '\t' << h.evaluate(s) << '\n';
    }
}

void testSSin(hermiteCubic h, int n, int pars) {
    h.d[0] = 1;
    h.d[pars - 1] = 0;
    smoothCubic_DEC(h);
    if(n == 0 || n == 4 || n == 8) {
        ofstream t3("Sine" + to_string(n + 1));
        for(int i = 0; i <= 3 * pars; i ++) {
            double s;
            s = i * M_PI_2 / (3 * pars);
            t3 << s << '\t' << h.evaluate(s) << '\n';
        }
        t3.close();
    }
}

void testSin(hermiteCubic h, int n, int pars) {
    not_a_knot(h);
    if (n == 0 || n == 4 || n == 8) {
        ofstream t2("Sine" + to_string(n));
        for (int i = 0; i <= 3 * pars; i++) {
            double s;
            s = i * M_PI_2 / (3 * pars);
            t2 << s << '\t' << h.evaluate(s) << '\n';
        }
        t2.close();
    }
}

int main() {
    //Modeling x^3 on [0,2] with one extra point of x=0.75
    /*uint n = 5;
    double dl = 0;
    double dr = 8.0;
    double p[5] = {0, 0.75, 1.25, 1.5, 2.0};
    double v[5] = {0, pow(0.75, 3), pow(1.25, 3), pow(1.5, 3), 8.0};
    double d[5];
    d[0] = dl;
    d[4] = dr;*/
    uint n[5] = {5, 10, 25, 100, 1000};
    double dl = 2;
    double dr = 10.4;
    double p[5] = {1, 1.5, 4, 4.6, 5.2};
    double v[5] = {1, pow(1.5,2), pow(4,2), pow(4.6, 2), pow(5.2, 2)};
    double d[5];
    d[0] = dl;
    d[4] = dr;
    hermiteCubic h(n[0], d, p, v);

    test1(h);
    //test2(h);

    double p2[7] = {1, 1.5, 3.2, 4, 4.6, 5.2, 5.5};
    double v2[7] = {1, pow(1.5,3), pow(3.2, 3), pow(4,3), pow(4.6, 3), pow(5.2, 3), pow(5.5, 3)};
    double d2[7];
    hermiteCubic h2(7, d2, p2, v2);
    //test1(h2);
    test2(h2);
    ofstream e1("SineErrorKnot");
    ofstream e2("SineErrorSmooth");
    //Sine test
    uint nsin[10] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    for(uint i = 0; i < 10; i ++) {
        double psin[2048];
        double vsin[2048];
        double dsin[2048];
        double dp = M_PI_2 / (nsin[i]);
        for(int j = 0; j < nsin[i]; j ++) {
            psin[j] = (dp * j);
            vsin[j] = (sin(psin[j]));
        }

        hermiteCubic hcsin(nsin[i], dsin, psin, vsin);
        testSSin(hcsin, i, nsin[i]);
        e2 << log10(nsin[i]) << '\t' << maxError(hcsin) << '\n';
        testSin(hcsin, i, nsin[i]);
        e1 << log10(nsin[i]) << '\t' << maxError(hcsin) << '\n';

    }
    e1.close();
    e2.close();

    system("gnuplot ./plot.gnu");
    //system("rm Sine* Part*");
    system("rm Sine0 Sine1 Sine4 Sine5 Sine8 Sine9 SineErrorSmooth SineErrorKnot Part1Curve Part2Curve");
    return 0;
}