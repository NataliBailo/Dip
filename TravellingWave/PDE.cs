using System;
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
        [Description("Количество точек по x")]
        public int N
        {
            get;
            set;
        }

        [Description("Количество точек по t")]
        public int M
        {
            get;
            set;
        }

        [Description("Интервал по x [-L, L]")]
        public double L
        {
            get;
            set;
        }

        [Description("Интервал по t [0, T]")]
        public double T
        {
            get;
            set;
        }

        [Description("Константа перед ядром")]
        public double b
        {
            get;
            set;
        }

        [Description("Константа перед U_xx")]
        public double D
        {
            get;
            set;
        }

        [Description("Задержка в Delta-Ядре")]
        public double d
        {
            get;
            set;
        }

        [Description("Константа 'a' в функции f(u) при Classical == false")]
        public double a
        {
            get;
            set;
        }

        [Description("Уравнение с классической функцией f(u)?")]
        public bool Classical
        {   // Уравнение с классической функцией f?
            get;
            set;
        }

        [Description("Delta-ядро?")]
        public bool DeltaCoupling
        {   // если да, то решается с дельта ядром
            get;
            set;
        }

        // Конструктор
        public PDE()
        {
            N = 2000;
            M = 2000;
            L = 50.0;
            T = 40.0;

            b = 1.0;
            d = 1.0;
            D = 1.0;
            a = 0.1;

            Classical = false;
            DeltaCoupling = true;
        }

        // методы
        public void load()
        {   // инициализируем/декларируем массивы и некоторые переменные
            hx = 2 * L / N; // шаг по x
            ht = T / M;  // шаг по t

            x = new double[N + 1];
            for (int i = 0; i < N + 1; i++) 
                x[i] = - L + i * hx;

            t = new double[M + 1];
            for (int j = 0; j < M + 1; j++) 
                t[j] = j * ht;

            u = new double[M + 1, N + 1];
        }

        public void initials()
        {   // Применяем начальные условия к функции u(x,t)

            for (int i = 0; i < N + 1; i++) // начальное значение при t = 0
                u[0, i] = u_x_0(x[i]);
        }

        public int solve()
        {   
            double step = ht / (hx * hx);

            double[] P = new double[N + 1];
            double[] Q = new double[N + 1];

            double[] ai = new double[3] { 0, -step * D, 1 };
            double[] bi = new double[3] { -1, -1 - 2 * step * D, 1 };
            double[] ci = new double[3] { -1, -step * D, 0 };
            double[] di = new double[N + 1];
            di[0] = 0; di[N] = 0; // если меняем условие Неймана (что-то вместо du/dn = 0), то эту строчку нужно закомментировать

            P[0] = ci[0] / bi[0];
            for (int i = 1; i < N; i++) 
                P[i] = ci[1] / (bi[1] - ai[1] * P[i - 1]);
            P[N] = ci[2] / (bi[2] - ai[2] * P[N - 1]);

            for (int j = 0; j < M; j++)
            {
                int extCode = progonkaJLayer(Q, P, ai, bi, di, j);
                if (extCode != 0) // если не получилось сделать прогонку, то прекращаем решать уравнение и показываем ошибку
                    return -1;
            }

            return 0;
        }

        private int progonkaJLayer(double[] Q, double[] P, double[] ai, double[] bi, double[] di, int j)
        {
            int k = Convert.ToInt32(d / hx);
            d = hx * k;

            //di[0] = ht * u_0_t(t[j]); // если условие Неймана ненулевое
            Q[0] = -di[0] / bi[0];
            for (int i = 1; i < N; i++)
            {
                calculateDCoeff(di, i, j, k);

                Q[i] = (ai[1] * Q[i - 1] - di[i]) / (bi[1] - ai[1] * P[i - 1]);

                // ловим, чтобы Q не равнялось бесконечности
                if (Double.IsNaN(Q[i]))
                    return -1;
            }
            //di[n] = ht * u_l_t(t[j]); // если условие Неймана ненулевое
            Q[N] = (ai[2] * Q[N - 1] - di[N]) / (bi[2] - ai[2] * P[N - 1]);

            u[j + 1, N] = Q[N];
            for (int i = N - 1; i > -1; i--) 
                u[j + 1, i] = P[i] * u[j + 1, i + 1] + Q[i];

            return 0;
        }

        private void calculateDCoeff(double[] di, int i, int j, int k)
        {
            if (DeltaCoupling)
            {
                if (b != 0)
                {
                    if (i - k <= 1) // если x - d <= -l
                    {
                        if (i + k <= N - 1) // если x - d <= -l И x + d <= l
                            di[i] = u[j, i] + ht * (b * (0.5 * (u[j, 1] + u[j, i + k]) - u[j, i]) + f(u[j, i]));
                        else if (i + k > N - 1) // если x - d <= -l И x + d > l
                            di[i] = u[j, i] + ht * (b * (0.5 * (u[j, 1] + u[j, N - 1]) - u[j, i]) + f(u[j, i]));
                    }
                    else if (i - k > 1) // если x - d > -l
                    {
                        if (i + k <= N - 1) // если x - d > -l И x + d <= l
                            di[i] = u[j, i] + ht * (b * (0.5 * (u[j, i - k] + u[j, i + k]) - u[j, i]) + f(u[j, i]));
                        else if (i + k > N - 1)  // если x - d > -l И x + d > l
                            di[i] = u[j, i] + ht * (b * (0.5 * (u[j, i - k] + u[j, N - 1]) - u[j, i]) + f(u[j, i]));
                    }
                }
                else
                    di[i] = u[j, i] + ht * (f(u[j, i]));
            }
            else if (b != 0)
                di[i] = u[j, i] + ht * (b * (integral(j, i) - u[j, i]) + f(u[j, i]));
            else
                di[i] = u[j, i] + ht * (f(u[j, i]));
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

        // разные функции
		private double f(double u)
        {
            if (Classical)
                return u - u * u * u / 3;
            else
                return -u * (u - 1) * (u - a);
        }

        private double integral(int j, int i)
        {   // Формула трапеции (для ядра K)
            // интегрируем в точке (t[j], x[i]) от -l до l
            double sum = 0;
            sum += kernel(x[i] - x[0]) * u[j, 0];
            for (int k = 1; k < N; k++) sum += 2 * kernel(x[i] - x[k]) * u[j, k];
            sum += kernel(x[i] - x[N]) * u[j, N];

            return hx * sum / 2;
        }

        private double kernel(double z)
        {
            //return Math.Exp(-z * z / 2) / Math.Sqrt(2 * Math.PI);
            return 0.5 * Math.Exp(-Math.Abs(z + 2));
        }            

        private double u_x_0(double x)
        {	// начальное условие при t = 0
            return 0.5 * (1 + Math.Tanh(x / (2 * Math.Sqrt(2))));
        }

        private double u_0_t(double t) { return 0.0; } // Условие Неймана на границе x = -l

        private double u_l_t(double t) { return 0.0; } // Условие Неймана на границе x = l
    }
}
