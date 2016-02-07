using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;

namespace TravellingWave
{
    class PDE
    {
        // переменные и массивы
        private double[] x, t;
        private double hx, ht; // шаги по x и по t
        private double[,] u;
        
        // свойства (properties)
        public int N
        {   // количество точек по x
            get;
            set;
        }

        public int M
        {   // количество точек по t
            get;
            set;
        }

        [Description("Interval for x [-L, L]")]
        public double L
        {   // bound x's segment
            get;
            set;
        }

        [Description("Interval for t [0, T]")]
        public double T
        {   // bound t's segment
            get;
            set;
        }

        [Description("constant in front of Kernel")]
        public double B
        {   // constant before coupling
            get;
            set;
        }

        [Description("Delay in Delta-Kernel")]
        public double D
        {   // a delay in delta Kernel
            get;
            set;
        }

        [Description("Current I excitatory")]
        public double I
        {   // current
            get;
            set;
        }

        [Description("F's constant in non-classical non-linearity (classical == false)")]
        public double A
        {   // f's constant if non-classical
            get;
            set;
        }

        [Description("Whether this is a classical f(u) or not")]
        public bool Classical
        {   // Is the equation with classical non-linearity?
            get;
            set;
        }

        [Description("Delta-Kernel or not?")]
        public bool DeltaCoupling
        {   // bool for deciding which equation solves
            get;
            set;
        }

        // Constructor with default values
        public PDE()
        {
            N = 2000;
            M = 2000;
            L = 50.0;
            T = 40.0;
            B = 0.0;
            D = 1.0;
            I = 0.0;
            A = 0.1;

            Classical = true;
            DeltaCoupling = true;
        }

        // methods
        public void load()
        {   // initialize/declare arrays and steps
            // If we want to change one of the parameters: n, m, l, TB,
            // then it needs to call this (plus Intiials) functions again.
            hx = 2 * L / N; // step for x
            ht = T / M;  // step for t

            x = new double[N + 1]; // arrange x's
            for (int i = 0; i < N + 1; i++) x[i] = - L + i * hx;

            t = new double[M + 1]; // arrange t's
            for (int j = 0; j < M + 1; j++) t[j] = j * ht;

            u = new double[M + 1, N + 1];
        }

        public void initials()
        {   // Initialize initials

            for (int i = 0; i < N + 1; i++) // setting an initial waves at t = 0
                u[0, i] = u_x_0(x[i]);
        }

        public int solve()
        {
            // If we changed ONLY eps, beta, gamma, b, d, Kernel, f or Iext,
            // then just recall this function.
            
            double step = ht / (hx * hx);

            double[] P = new double[N + 1];
            double[] Q = new double[N + 1];

            double[] ai = new double[3] { 0, -step, 1 };
            double[] bi = new double[3] { -1, -1 - 2 * step, 1 };
            double[] ci = new double[3] { -1, -step, 0 };
            double[] di = new double[N + 1];
            di[0] = 0; di[N] = 0; // if Neumann condition changes (smth except du/dn = zero), it needs to be commented

            P[0] = ci[0] / bi[0];
            for (int i = 1; i < N; i++) P[i] = ci[1] / (bi[1] - ai[1] * P[i - 1]);
            P[N] = ci[2] / (bi[2] - ai[2] * P[N - 1]);

            for (int j = 0; j < M; j++)
            {
                int extCode = progonkaJLayer(Q, P, ai, bi, di, j);
                if (extCode != 0)
                    return -1;
            }

            return 0;
        }

        private int progonkaJLayer(double[] Q, double[] P, double[] ai, double[] bi, double[] di, int j)
        {
            int k = Convert.ToInt32(D / hx);
            D = hx * k;

            //di[0] = ht * u_0_t(t[j]); // if Neumann condition is not a zero
            Q[0] = -di[0] / bi[0];
            for (int i = 1; i < N; i++)
            {
                calculateDCoeff(di, i, j, k);

                Q[i] = (ai[1] * Q[i - 1] - di[i]) / (bi[1] - ai[1] * P[i - 1]);

                // catching Q is NaN
                if (Double.IsNaN(Q[i]))
                    return -1;
            }
            //di[n] = ht * u_l_t(t[j]); // if Neumann condition is not a zero
            Q[N] = (ai[2] * Q[N - 1] - di[N]) / (bi[2] - ai[2] * P[N - 1]);

            u[j + 1, N] = Q[N];
            for (int i = N - 1; i > -1; i--) u[j + 1, i] = P[i] * u[j + 1, i + 1] + Q[i];

            return 0;
        }

        private void calculateDCoeff(double[] di, int i, int j, int k)
        {
            if (DeltaCoupling)
            {
                if (B != 0)
                {
                    if (i - k <= 1) // if x - d <= -l
                    {
                        if (i + k <= N - 1) // if x - d <= -l and x + d <= l
                            di[i] = u[j, i] + ht * (B * (u[j, 1] + u[j, i + k] - 2 * u[j, i]) + f(u[j, i]) + I);
                        else if (i + k > N - 1) // if x - d <= -l and x + d > l
                            di[i] = u[j, i] + ht * (B * (u[j, 1] + u[j, N - 1] - 2 * u[j, i]) + f(u[j, i]) + I);
                    }
                    else if (i - k > 1) // if x - d > -l
                    {
                        if (i + k <= N - 1) // if x - d > -l and x + d <= l
                            di[i] = u[j, i] + ht * (B * (u[j, i - k] + u[j, i + k] - 2 * u[j, i]) + f(u[j, i]) + I);
                        else if (i + k > N - 1)  // if x - d > -l and x + d > l
                            di[i] = u[j, i] + ht * (B * (u[j, i - k] + u[j, N - 1] - 2 * u[j, i]) + f(u[j, i]) + I);
                    }
                }
                else
                    di[i] = u[j, i] + ht * (f(u[j, i]) + I);
            }
            else if (B != 0)
                di[i] = u[j, i] + ht * (B * (integral(j, i) - u[j, i]) + f(u[j, i]) + I);
            else
                di[i] = u[j, i] + ht * (f(u[j, i]) + I);
        }

        public double getX(int i)
        { 
            return x[i]; 
        }

        public double getT(int j)
        {
            return t[j];
        }
        
        public double getU(int j, int i)
        { 
            return u[j, i]; 
        }

        // various functions
		private double f(double u)
        {
            if (Classical)
                return u - u * u * u / 3;
            else
                return -u * (u - 1) * (u - A);
        }

        private double integral(int j, int i)
        {   // Trapezoidal rule, uniform grid
            // integrating at point (t[j], x[i]) from -l to l
            double sum = 0;
            sum += kernel(x[i] - x[0]) * u[j, 0];
            for (int k = 1; k < N; k++) sum += 2 * kernel(x[i] - x[k]) * u[j, k];
            sum += kernel(x[i] - x[N]) * u[j, N];

            return hx * sum / 2;
        }

        private double kernel(double z)
        {
            //return Math.Exp(-z * z / 2) / Math.Sqrt(2 * Math.PI);
            return 1.0 / 2 * Math.Exp(-Math.Abs(z + 2));
        }            

        private double u_x_0(double x)
        {	// начальное условие при t = 0
            return (1.0 / 2) * (1 + Math.Tanh(x / (2 * Math.Sqrt(2))));
        }

        private double u_0_t(double t) { return 0.0; } // Neumann boundary condition at x = -l

        private double u_l_t(double t) { return 0.0; } // Neumann boundary condition at x = l
    }
}
