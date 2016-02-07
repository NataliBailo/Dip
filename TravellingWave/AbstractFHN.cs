using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;

namespace TravellingWave
{
    abstract class AbstractFHN
    {
        protected double hx, ht; // steps
        protected double[] t; // time

        public AbstractFHN() {
            Eps = 0.08;
            Beta = 0.064;
            Alpha = 0.056;
            A = 0.1;
            Classical = true;
        }

        public abstract void load();
        public abstract void initials();
        public abstract int solve();

        public int N
        {   // quantity of u,v t's (x's)
            get;
            set;
        }

		public virtual double L
        {
            get;
            set;
        }

        [Description("Interval for t [0, T]")]
        public double T
        {   // bound t's segment
            get;
            set;
        }

        [Description("v's equation constant")]
        public double Eps
        {   // v's equation constants
            get;
            set;
        }

        [Description("v's equation constant")]
        public double Alpha
        {   // v's equation constants
            get;
            set;
        }

        [Description("v's equation constant")]
        public double Beta
        {   // v's equation constants
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

        public double getT(int j)
        {
            return t[j];
        }

        protected double f(double u)
        {
            if (Classical)
                return u - u * u * u / 3;
            else
                return -u * (u - 1) * (u - A);
        }
    }
}
