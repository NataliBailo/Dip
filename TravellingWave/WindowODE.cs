using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace TravellingWave
{
    public partial class WindowODE : Form
    {
        ODE[] odes;

        public WindowODE()
        {
            InitializeComponent();
        }

        private void WindowODE_Load(object sender, EventArgs e)
        {
            loadEquations(1);
        }

        private void loadEquations(int num)
        {
            odes = new ODE[num];

            for (int i = 0; i < num; i++)
                odes[i] = new ODE();

            propertyGrid1.SelectedObject = odes[0];

            if (num == 2)
                propertyGrid2.SelectedObject = odes[1];
            else
                propertyGrid2.SelectedObject = null;
        }

        private void plot(int j, ODE obj, int numEq)
        {   // plots single point on t plot and phase plane
            chart.Series[2 * numEq].Points.AddXY(obj.getT(j), obj.getU(j));
            chart.Series[2 * numEq + 1].Points.AddXY(obj.getT(j), obj.getV(j));

            chartPhase.Series[3 * numEq].Points.AddXY(obj.getU(j), obj.getV(j));
        }

        private void plotNullclines(ODE obj, int numEq)
        {
            for (int j = 0; j < obj.N + 1; j++)
            {   // drawing nullclines
                chartPhase.Series[3 * numEq + 1].Points.AddXY(obj.getUN(j), obj.getV1(j));
                chartPhase.Series[3 * numEq + 2].Points.AddXY(obj.getUN(j), obj.getV2(j));
            }
        }

        private void propertyGrid1_SelectedGridItemChanged(object sender, SelectedGridItemChangedEventArgs e)
        {
            btnSolveBeh();
        }

        private void propertyGrid2_SelectedGridItemChanged(object sender, SelectedGridItemChangedEventArgs e)
        {
            btnSolveBeh();
        }

        private void btnSolveBeh()
        {
            // Solve button behaviour
            // if we change I, tau, a, b, Kernel (?), f or initials,
            // then disable Plot and call Solve again.
            if (btnPlot.Enabled)
            {
                btnPlot.Enabled = false;
                btnSolve.Enabled = true;
            }

            lblError.Visible = false;

            trBarT.Value = 0;
            trBarT.Enabled = false;
            timerT.Enabled = false;
        }

        private void btnPlotBeh()
        {
            // Plot button behaviour
            if (!lblError.Visible)
            {
                btnSolve.Enabled = false;
                btnPlot.Enabled = true;

                trBarT.Value = 0;
                trBarT.Enabled = true;
            }
        }

        private void btnSolve_Click(object sender, EventArgs e)
        {
            prBarSolve.Value = 0;
            prBarSolve.Maximum = 3;
            trBarT.Maximum = odes[0].N;

            for (int i = 0; i < odes.Length; i++)
                odes[i].load();
            prBarSolve.Value++;

            for (int i = 0; i < odes.Length; i++)
                odes[i].initials();
            prBarSolve.Value++;

            for (int i = 0; i < odes.Length; i++)
            {
                int extCode = odes[i].solve();
                if (extCode != 0)
                    lblError.Visible = true;
            }
            prBarSolve.Value++;

            btnPlotBeh();
        }

        private void btnPlot_Click(object sender, EventArgs e)
        {
            setPlot();

            // if 'Plot all' checked then need to clear plot.
            if (rdBtnPlotAll.Checked)
            {
                //trBarT.Enabled = false;
                trBarT.Value = 0;

                clearPlot();

                for (int i = 0; i < odes.Length; i++)
                    plotNullclines(odes[i], i);

                for (int i = 0; i < odes.Length; i++)
                {
                    for (int j = 0; j < odes[i].N + 1; j++)
                    {
                        plot(j, odes[i], i);
                    }
                }
            }

            if (rdBtnTmr.Checked)
            {
                clearAllButNullclinesPlot();

                trBarT.Enabled = true;
                timerT.Enabled = true;
            }
            else
            {
                //trBarT.Enabled = false; 
                timerT.Enabled = false;
            }
        }

        private void timerT_Tick(object sender, EventArgs e)
        {
            if (trBarT.Value == odes[0].N)
            {
                for(int i =0 ; i < odes.Length; i++)
                    plot(trBarT.Value, odes[i], i);

                trBarT.Value = 0;

                clearAllButNullclinesPlot();
            }
            else
            {
                trBarT.Value++;
                for (int i = 0; i < odes.Length; i++)
                    plot(trBarT.Value, odes[i], i);
            }
        }

        private void trBarT_Scroll(object sender, EventArgs e)
        {
            clearAllButNullclinesPlot();
            timerT.Enabled = false;

            // plot all 'till moment on the trackbar
            for (int j = 0; j < trBarT.Value; j++)
            {
                for (int i = 0; i < odes.Length; i++)
                    plot(j, odes[i], i);
            }
        }

        private void btnTuneT_Click(object sender, EventArgs e)
        {
            chart.ChartAreas[0].AxisY.Maximum = Convert.ToDouble(txtBoxMaxUVT.Text);
            chart.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUVT.Text);
        }

        private void btnTunePhase_Click(object sender, EventArgs e)
        {
            chartPhase.ChartAreas[0].AxisY.Maximum = Convert.ToDouble(txtBoxMaxUVPhase.Text);
            chartPhase.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUVPhase.Text);
        }

        private void btnStopTimer_Click(object sender, EventArgs e)
        {
            if (timerT.Enabled)
            {
                rdBtnTmr.Checked = false;
                timerT.Enabled = false;
                //trBarT.Enabled = false;

                rdBtnPlotAll.Checked = true;
            }
        }

        private void clearAllButNullclinesPlot()
        {
            for (int i = 0; i < chart.Series.Count(); i++)
                chart.Series[i].Points.Clear();

            for (int i = 0; i < 2; i++)
                chartPhase.Series[3 * i].Points.Clear();
        }

        private void clearPlot()
        {
            for (int i = 0; i < chart.Series.Count(); i++)
                chart.Series[i].Points.Clear();

            for (int i = 0; i < chartPhase.Series.Count(); i++)
                chartPhase.Series[i].Points.Clear();
        }

        private void setPlot()
        {
            chart.ChartAreas[0].AxisX.Minimum = 0;
            chart.ChartAreas[0].AxisX.Maximum = odes[0].T + 1;

            chart.ChartAreas[0].AxisY.Maximum = Convert.ToDouble(txtBoxMaxUVT.Text);
            chart.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUVT.Text);

            chart.ChartAreas[0].AxisX.Interval = Convert.ToInt32((chart.ChartAreas[0].AxisX.Maximum + chart.ChartAreas[0].AxisX.Minimum) / 5.0);
            chart.ChartAreas[0].AxisY.Interval = Convert.ToInt32((chart.ChartAreas[0].AxisY.Maximum + chart.ChartAreas[0].AxisY.Minimum) / 4.0);

            chartPhase.ChartAreas[0].AxisX.Minimum = -odes[0].L;
            chartPhase.ChartAreas[0].AxisX.Maximum = odes[0].L;

            chartPhase.ChartAreas[0].AxisY.Maximum = 1.5 + odes[0].I;
            chartPhase.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUVPhase.Text);

            chart.Series[2].Color = Color.Blue;
            chart.Series[3].Color = Color.OrangeRed;

            chartPhase.Series[3].Color = Color.Blue;
        }

        private void checkBox2ndEq_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox2ndEq.Checked)
                loadEquations(2);
            else
                loadEquations(1);
        }

        private void btnAbout_Click(object sender, EventArgs e)
        {
            AboutODE o = new AboutODE();
            o.Show();
        }

        private void btnTickSlower_Click(object sender, EventArgs e)
        {
            timerT.Interval += odes[0].N / 5;
        }

        private void btnTickFaster_Click(object sender, EventArgs e)
        {
            if (timerT.Interval - odes[0].N / 10 > 0)
                timerT.Interval -= odes[0].N / 10;
            else
            {
                timerT.Interval = 1;
            }
        }
    }
}
