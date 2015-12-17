// #define forn(i,n) for (int i=0;i<n;i++)
// #define rforn(i,n) for (int i=n-1;i>=0;i--)
// #define mp(a,b) make_pair(a,b)
// #define LL long long
// #define N 400002
// #define MOD 1000000007
#include <iostream>
#include <math.h>
#include <stdio.h>
#define VNAME(x) #x

using namespace std;
inline void create(double **mtr, int n)
{
    for(int i=0; i<n; ++i)
    {
        mtr[i] = new double[n]; memset(mtr[i], 0, sizeof(double) * n);
    }
}
inline void print(double **mtr, int n)
{
    FILE * file = fopen("/Users/ramil/workspace/study/lapin_tasks/output_for_matrix.txt","w");
    for(int i=0; i<n; ++i)
    {
        fprintf(file, "matrix name: %s\n", VNAME(mtr));
        for(int j=0; j<n; ++j)
          fprintf(file, "%.6f ", mtr[i][j]);
        fprintf(file, "\n");
    }
}
inline double f(int i, int j, double h)
{
    double x = i*h - 0.5, y = j*h - 0.5; if(x*x + y*y <= 0.1)
    return 0; return 4;
}
int solve(double h, double sigma, double epsilon, double tao)
{
    FILE * file = fopen("/Users/ramil/workspace/study/lapin_tasks/output.txt","w");

    int m = round(1/h);
    int n = m+1;
    double h2 = h*h;
    double **u = new double*[n], **F = new double*[n], **p1 = new double*[n], **p2 = new double*[n], **l1 = new double*[n], **l2 = new double*[n];
    create(u, n); create(F, n); create(p1, n); create(p2, n); create(l1, n); create(l2, n); int iter=0;
    for(;; ++iter)
    {
        for(int i=0; i<n; ++i)
        {
            for(int j=0; j<n; ++j)
            {
                double modl = sqrt(l1[i][j] * l1[i][j] + l2[i][j] * l2[i][j]); if(modl <= 1) p1[i][j] = l1[i][j], p2[i][j] = l2[i][j];
                else p1[i][j] = l1[i][j]/modl, p2[i][j] = l2[i][j]/modl;
            }
        }
        for(int i=1; i<m; ++i)
        {
            for(int j=1; j<m; ++j)
            {
                double Su = (2 * (u[i][j] - u[i-1][j]) + 0 * (u[i][j] - u[i][j-1])) / h;
                double Lp = 2 * (p1[i][j] - p1[i+1][j] + p2[i][j] - p2[i][j+1]) / h - (l1[i][j] - l1[i+1][j] + l2[i][j] - l2[i][j+1]) / h;
                F[i][j] = Su/2 + Lp/2 + f(i,j,h) / 2;
            }
        }
        for(int it=0;;++it)
        {
            for(int i=1; i<m; ++i)
            for(int j=1; j<m; ++j)
            u[i][j] = (1 - sigma) * u[i][j] + sigma * 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] + h2 *
            F[i][j]);
            double norm_r = 0;
            for(int i=1; i<m; ++i)
            {
                for(int j=1; j<m; ++j)
                {
                    double r = (4 * u[i][j] - u[i-1][j] - u[i+1][j] - u[i][j-1] - u[i][j+1]) / h2 - F[i][j];
                    norm_r += h2 * r * r;
                }
            }
            if(sqrt(norm_r) < epsilon)
            {
                fprintf(file, "%s = %d\n",VNAME(it), it);
                break;
            }
        }
        for(int i=1; i<n; ++i)
        {
            for(int j=1; j<n; ++j)
            {
                l1[i][j] = l1[i][j] - tao * (p1[i][j] - (u[i][j] - u[i-1][j]) / h); l2[i][j] = l2[i][j] - tao * (p2[i][j] - (u[i][j] - u[i][j-1]) / h);
            }
        }
        double r = 0;
        for(int i=1; i<n; ++i)
        {
            for(int j=1; j<n; ++j)
            {
                double d1 = p1[i][j] - (u[i][j] - u[i-1][j]) / h,
                d2 = p2[i][j] - (u[i][j] - u[i][j-1]) / h; r += d1 * d1 + d2 * d2;
            }
        }
        r = h * sqrt(r);
        fprintf(file, "%s %.5f\n",VNAME(r), r);
        if(r < epsilon) break;
    }
    fprintf(file, "\n");
    print(u, n);
    for(int i=0; i<n; ++i)
    {
        fprintf(file, "\n");
        for(int j=0; j<n; ++j)
        {
            fprintf(file, "sqrt of squares = %.6f ", sqrt(p1[i][j] * p1[i][j] + p2[i][j] * p2[i][j]));
        }
    }
    return iter;
}
int main()
{
    solve(0.01, 1.92, 0.01, 1.9);
    return 0;
}
