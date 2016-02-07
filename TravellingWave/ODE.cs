using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;

namespace TravellingWave
{
    class ODE : AbstractFHN
    {
        // variables and arrays
        private double[] u, v;
        private double[] u_null, v1, v2; // nullclines

        // properties
        public double U0
        {   // intiial u
            get;
            set;
        }

        public double V0
        {   // intiial v
            get;
            set;
        }

        // Constructor with default parameters
        public ODE() : base()
        {
            N = 1000;
            L = 2.5;
            T = 100.0;
            U0 = 1.0;
            V0 = 0.1;
            I = 0.5;
        }

        // methods
        public override void load()
        {   // initialize/declare arrays and steps
            // If we want to change one of the parameters: n or TB,
            // then it needs to call this (plus Intiials) functions again.

            hx = 2 * L / N; //
            ht = T / N;  // step for t

            t = new double[N + 1]; // arrange t's
            for (int j = 0; j < N + 1; j++) 
                t[j] = j * ht;

            u = new double[N + 1];
            v = new double[N + 1];

            v1 = new double[N + 1];
            v2 = new double[N + 1];
            
            u_null = new double[N + 1];
            for (int j = 0; j < N + 1; j++) 
                u_null[j] = - L + j * hx;
        }

        public override void initials()
        {   // Initialize initials
            u[0] = U0;
            v[0] = V0;
        }

        public override int solve()
        {   // If we changed ONLY alpha, beta, Iext, Kernel or f (either a),
            // then just recall this function.

            double utemp, vtemp;

            for (int j = 0; j < N; j++)
            {
                utemp = u[j] + ht * f1(u[j], v[j]);
                vtemp = v[j] + ht * f2(u[j], v[j]);

                if (Double.IsNaN(utemp) || Double.IsNaN(vtemp))
                    return -1;

                //// Euler's 1st order
                //u[j + 1] = utemp;
                //v[j + 1] = vtemp;

                // Euler's 2nd order
                u[j + 1] = u[j] + ht / 2 * (f1(u[j], v[j]) + f1(utemp, vtemp));
                v[j + 1] = v[j] + ht / 2 * (f2(u[j], v[j]) + f2(utemp, vtemp));
            }
            
            nullclines();

            return 0;
        }

        public void nullclines()
        {
            if (Beta != 0.0)
            {
                for (int j = 0; j < N + 1; j++)
                {
                    v1[j] = f(u_null[j]) + I;
                    v2[j] = (u_null[j] + Eps) / Beta;
                }
            }
        }

        public double getU(int j)
        { 
            return u[j]; 
        }

        public double getV(int j)
        { 
            return v[j]; 
        }

        public double getUN(int j)
        { 
            return u_null[j]; 
        }

        public double getV1(int j)
        { 
            return v1[j]; 
        }

        public double getV2(int j)
        { 
            return v2[j]; 
        }

        // various functions
        private double f1(double u, double v)
        {
            return f(u) - v + I;
        }

        private double f2(double u, double v)
        {
            return Eps * u + Alpha - Beta * v;
        }
    }
}
