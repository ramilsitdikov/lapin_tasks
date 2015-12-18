#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
#define VNAME(x) #x
#define PI 3.14159265

using namespace std;
inline void create(double **mtr, int n){
    for(int i=0; i<n; ++i){ mtr[i] = new double[n]; memset(mtr[i], 0, sizeof(double) * n); }
}

inline void print_matrix(double **mtr, int n, string path){
    FILE * file = fopen(path.c_str(), "w");
    fprintf(file, "[");
    for(int i=0; i<n; ++i){ fprintf(file, "[");
        for(int j=0; j<n; ++j){ fprintf(file, "%.6f", mtr[i][j]); if (j != n - 1) fprintf(file, ","); }
            fprintf(file, "]");
            if (i != n - 1) fprintf(file, ",");}
    fprintf(file, "]");
}

inline void print_list(vector <float> mtr, int n, string path){
    FILE * file = fopen(path.c_str(),"w");
    fprintf(file, "[");
    for(int j=0; j<n; ++j){
        fprintf(file, "%.6f", mtr[j]);
        if (j != n - 1) fprintf(file, ","); }
    fprintf(file, "]");
}

inline double f(int i, int j, double h){
    double x = i*h - 0.5, y = j*h - 0.5;
    if(x*x + y*y <= 0.1) return 1; return 4;
}

int solve(double h, double sigma, double epsilon, double tao, double k1, double k2){
    FILE * file = fopen("/Users/ramil/workspace/study/lapin_tasks/output.txt","w");

    int m = round(1/h);
    int n = m+1;
    double h2 = h*h;
    double **u = new double*[n],
                **F = new double*[n],
                **p1 = new double*[n],
                **p2 = new double*[n],
                **l1 = new double*[n],
                **l2 = new double*[n],
                **p = new double*[n];
    create(u, n); create(F, n); create(p1, n); create(p2, n); create(l1, n); create(l2, n); create(p, n);

    //start of big loop. Out if r < epsilon
    int iter=0;
    vector <float> iter_list;
    vector <float> r_list;
    for(;; ++iter){
        iter_list.resize(iter + 1);
        r_list.resize(iter + 1);

        for(int i=0; i<n; ++i){ for(int j=0; j<n; ++j){
                double modl = sqrt(l1[i][j] * l1[i][j] + l2[i][j] * l2[i][j]);
                if(modl <= 1) p1[i][j] = l1[i][j], p2[i][j] = l2[i][j];
                else p1[i][j] = l1[i][j]/modl, p2[i][j] = l2[i][j]/modl; } }
        for(int i=1; i<m; ++i){ for(int j=1; j<m; ++j){
                double Su = (k1 * (u[i][j] - u[i-1][j]) + k2 * (u[i][j] - u[i][j-1])) / h;
                double Lp = 2 * (p1[i][j] - p1[i+1][j] + p2[i][j] - p2[i][j+1]) / h - (l1[i][j] - l1[i+1][j] + l2[i][j] - l2[i][j+1]) / h;
                F[i][j] = Su/2 + Lp/2 + f(i,j,h) / 2; } }
        //Start of small loop. Out if sqrt(norm_r) < epsilon
        int it=0;
        for(;;++it) {
            for(int i=1; i<m; ++i) for(int j=1; j<m; ++j) u[i][j] = (1 - sigma) * u[i][j] + sigma * 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] + h2 * F[i][j]);
            double norm_r = 0;
            for(int i=1; i<m; ++i){ for(int j=1; j<m; ++j){
                    double r = (4 * u[i][j] - u[i-1][j] - u[i+1][j] - u[i][j-1] - u[i][j+1]) / h2 - F[i][j];
                    norm_r += h2 * r * r; } }
            if(sqrt(norm_r) < epsilon){ break; } }
        //end of small loop with it
        for(int i=1; i<n; ++i){ for(int j=1; j<n; ++j){
                l1[i][j] = l1[i][j] - tao * (p1[i][j] - (u[i][j] - u[i-1][j]) / h);
                l2[i][j] = l2[i][j] - tao * (p2[i][j] - (u[i][j] - u[i][j-1]) / h); } }
        double r = 0;
        for(int i=1; i<n; ++i){ for(int j=1; j<n; ++j){
                double d1 = p1[i][j] - (u[i][j] - u[i-1][j]) / h,
                            d2 = p2[i][j] - (u[i][j] - u[i][j-1]) / h;
                r += d1 * d1 + d2 * d2; } }
        r = h * sqrt(r);
        r_list[iter] = r;
        iter_list[iter] = it;
        fprintf(file, "Number of iteration: %d, amount of iterations = %d, residual norm is equal to: %f\n",iter + 1, it + 1, r);
        if(r < epsilon) break; }

    for(int i=0; i<n; ++i) for(int j=0; j<n; ++j) p[i][j] = sqrt(p1[i][j] * p1[i][j] + p2[i][j] * p2[i][j]);
    //end of big loop with iter

    //print out results
    string path_for_u = "/Users/ramil/workspace/study/lapin_tasks/u_matrix.txt";
    string path_for_p = "/Users/ramil/workspace/study/lapin_tasks/p_matrix.txt";
    string path_for_r = "/Users/ramil/workspace/study/lapin_tasks/r_list.txt";
    string path_for_iter = "/Users/ramil/workspace/study/lapin_tasks/iter_list.txt";

    fprintf(file, "\nNumber of iteration at all = %d\n", iter + 1);

    print_matrix(u, n, path_for_u);
    print_matrix(p, n, path_for_p);
    print_list(r_list, iter, path_for_r);
    print_list(iter_list, iter, path_for_iter);

    return 0;
}

int main()
{
    // int solve(double h, double sigma, double epsilon, double tao, double k1, double k2)
    solve(0.01, 1.92, 0.01, 2, 0, 2);
    return 0;
}
