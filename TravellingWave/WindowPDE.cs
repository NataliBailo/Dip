﻿using System;
using System.Linq;
using System.Windows.Forms;

namespace TravellingWave
{
    public partial class WindowPDE : Form
    {
        PDE[] pdes;

        public WindowPDE()
        {
            InitializeComponent();
        }

        private void WindowPDE_Load(object sender, EventArgs e)
        {
            loadEquations(1);
        }

        private void loadEquations(int num)
        {
            pdes = new PDE[num];

            for (int i = 0; i < num; i++)
                pdes[i] = new PDE();

            propertyGrid1.SelectedObject = pdes[0];

            if (num == 2)
                propertyGrid2.SelectedObject = pdes[1];
            else
                propertyGrid2.SelectedObject = null;
        }

        private void plot(int j, PDE obj, int numEq)
        {   // рисует весь слой по x при фиксированном t
            for (int i = 0; i < obj.N + 1; i++)
                chart.Series[numEq].Points.AddXY(obj.getX(i), obj.getU(j, i));
        }

        private void propertyGrid_PropertyValueChanged(object s, PropertyValueChangedEventArgs e)
        {
            btnSolveBeh();
        }

        private void btnSolveBeh()
        {   // вызвать, если меняем любой из параметров уравнения

            // выключаем кнопку "Нарисовать"
            if (btnPlot.Enabled)
            {
                btnPlot.Enabled = false;
                btnSolve.Enabled = true;
            }

            lblError.Visible = false;

            // сбрасываем значение трекбара
            // и выключаем таймер
            trBarT.Value = 0;
            trBarT.Enabled = false;
            timerT.Enabled = false;
        }

        private void btnPlotBeh()
        {
            // если не появилось ошибки,
            // то включаем кнопку "Нарисовать" и трекбар
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
            trBarT.Maximum = pdes[0].M;

            for (int i = 0; i < pdes.Length; i++)
                pdes[i].load();
            prBarSolve.Value++;

            for (int i = 0; i < pdes.Length; i++)
                pdes[i].initials();
            prBarSolve.Value++;

            for (int i = 0; i < pdes.Length; i++)
            {
                int extCode = pdes[i].solve();
                if (extCode != 0)
                    lblError.Visible = true;
            }
            prBarSolve.Value++;

            btnPlotBeh();
        }

        private void btnPlot_Click(object sender, EventArgs e)
        {
            clearPlot();
            setPlot();

            for (int i = 0; i < pdes.Length; i++)
                plot(trBarT.Value, pdes[i], i);

            if (rdBtnTmr.Checked)
                timerT.Enabled = true;
            else
                timerT.Enabled = false;
        }

        private void timerT_Tick(object sender, EventArgs e)
        {
            clearPlot();

            if (trBarT.Value == pdes[0].M)
            {   // если это последний сегмент по t, 
                // то рисуем его и сбрасываем значение трекбара на 0
                for (int i = 0; i < pdes.Length; i++)
                    plot(trBarT.Value, pdes[i], i);

                trBarT.Value = 0;
            }
            else
            {
                trBarT.Value++;
                for (int i = 0; i < pdes.Length; i++)
                    plot(trBarT.Value, pdes[i], i);
            }
        }

        private void trBarT_Scroll(object sender, EventArgs e)
        {
            clearPlot();

            for (int i = 0; i < pdes.Length; i++)
                plot(trBarT.Value, pdes[i], i);
        }

        private void btnTune_Click(object sender, EventArgs e)
        {
            try
            {
                chart.ChartAreas[0].AxisY.Maximum = Convert.ToDouble(txtBoxMaxUV.Text);
                chart.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUV.Text);
            }
            catch (Exception)
            {
                MessageBox.Show(
                    "Введено недопустимое значение (в качестве разделителя используйте запятую)",
                    "Недопустимое значение",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Exclamation);
            }
        }

        private void btnStopTimer_Click(object sender, EventArgs e)
        {
            if (timerT.Enabled)
            {
                rdBtnTmr.Checked = false;
                timerT.Enabled = false;
            }
        }

        private void clearPlot()
        {
            for (int i = 0; i < chart.Series.Count(); i++)
                chart.Series[i].Points.Clear();
        }

        private void setPlot()
        {
            chart.ChartAreas[0].AxisX.Minimum = -Convert.ToDouble(pdes[0].L);
            chart.ChartAreas[0].AxisX.Maximum = Convert.ToDouble(pdes[0].L);

            chart.ChartAreas[0].AxisY.Maximum = Convert.ToDouble(txtBoxMaxUV.Text);
            chart.ChartAreas[0].AxisY.Minimum = Convert.ToDouble(txtBoxMinUV.Text);

            chart.ChartAreas[0].AxisX.Interval = Convert.ToInt32((chart.ChartAreas[0].AxisX.Maximum + chart.ChartAreas[0].AxisX.Minimum) / 6.0);
            chart.ChartAreas[0].AxisY.Interval = Convert.ToInt32((chart.ChartAreas[0].AxisY.Maximum + chart.ChartAreas[0].AxisY.Minimum) / 6.0);
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
            AboutPDE o = new AboutPDE();
            o.Show();
        }
    }
}
