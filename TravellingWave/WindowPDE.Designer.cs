﻿namespace TravellingWave
{
    partial class WindowPDE
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea6 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend6 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series11 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series12 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.lblError = new System.Windows.Forms.Label();
            this.checkBox2ndEq = new System.Windows.Forms.CheckBox();
            this.lblMinUV = new System.Windows.Forms.Label();
            this.lblMaxUV = new System.Windows.Forms.Label();
            this.btnTune = new System.Windows.Forms.Button();
            this.txtBoxMinUV = new System.Windows.Forms.TextBox();
            this.txtBoxMaxUV = new System.Windows.Forms.TextBox();
            this.chart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.btnSolve = new System.Windows.Forms.Button();
            this.btnPlot = new System.Windows.Forms.Button();
            this.trBarT = new System.Windows.Forms.TrackBar();
            this.prBarSolve = new System.Windows.Forms.ProgressBar();
            this.rdBtnTmr = new System.Windows.Forms.RadioButton();
            this.btnStopTimer = new System.Windows.Forms.Button();
            this.timerT = new System.Windows.Forms.Timer(this.components);
            this.btnAbout = new System.Windows.Forms.Button();
            this.propertyGrid1 = new System.Windows.Forms.PropertyGrid();
            this.tabControl1 = new System.Windows.Forms.TabControl();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.tabPage2 = new System.Windows.Forms.TabPage();
            this.propertyGrid2 = new System.Windows.Forms.PropertyGrid();
            ((System.ComponentModel.ISupportInitialize)(this.chart)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.trBarT)).BeginInit();
            this.tabControl1.SuspendLayout();
            this.tabPage1.SuspendLayout();
            this.tabPage2.SuspendLayout();
            this.SuspendLayout();
            // 
            // lblError
            // 
            this.lblError.AutoSize = true;
            this.lblError.Font = new System.Drawing.Font("Microsoft Sans Serif", 9.75F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.lblError.ForeColor = System.Drawing.Color.DarkRed;
            this.lblError.Location = new System.Drawing.Point(21, 686);
            this.lblError.Name = "lblError";
            this.lblError.Size = new System.Drawing.Size(158, 16);
            this.lblError.TabIndex = 80;
            this.lblError.Text = "Ошибка вычислений";
            this.lblError.Visible = false;
            // 
            // checkBox2ndEq
            // 
            this.checkBox2ndEq.AutoSize = true;
            this.checkBox2ndEq.Location = new System.Drawing.Point(16, 442);
            this.checkBox2ndEq.Name = "checkBox2ndEq";
            this.checkBox2ndEq.Size = new System.Drawing.Size(141, 17);
            this.checkBox2ndEq.TabIndex = 82;
            this.checkBox2ndEq.Text = "Решать два уравнения";
            this.checkBox2ndEq.UseVisualStyleBackColor = true;
            this.checkBox2ndEq.CheckedChanged += new System.EventHandler(this.checkBox2ndEq_CheckedChanged);
            // 
            // lblMinUV
            // 
            this.lblMinUV.AutoSize = true;
            this.lblMinUV.Location = new System.Drawing.Point(202, 607);
            this.lblMinUV.Name = "lblMinUV";
            this.lblMinUV.Size = new System.Drawing.Size(41, 13);
            this.lblMinUV.TabIndex = 67;
            this.lblMinUV.Text = "min u,v";
            // 
            // lblMaxUV
            // 
            this.lblMaxUV.AutoSize = true;
            this.lblMaxUV.Location = new System.Drawing.Point(199, 40);
            this.lblMaxUV.Name = "lblMaxUV";
            this.lblMaxUV.Size = new System.Drawing.Size(44, 13);
            this.lblMaxUV.TabIndex = 66;
            this.lblMaxUV.Text = "max u,v";
            // 
            // btnTune
            // 
            this.btnTune.Location = new System.Drawing.Point(242, 333);
            this.btnTune.Name = "btnTune";
            this.btnTune.Size = new System.Drawing.Size(42, 23);
            this.btnTune.TabIndex = 73;
            this.btnTune.Text = "Zoom";
            this.btnTune.UseVisualStyleBackColor = true;
            this.btnTune.Click += new System.EventHandler(this.btnTune_Click);
            // 
            // txtBoxMinUV
            // 
            this.txtBoxMinUV.Location = new System.Drawing.Point(243, 604);
            this.txtBoxMinUV.Name = "txtBoxMinUV";
            this.txtBoxMinUV.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMinUV.TabIndex = 57;
            this.txtBoxMinUV.Text = "-1,5";
            this.txtBoxMinUV.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // txtBoxMaxUV
            // 
            this.txtBoxMaxUV.Location = new System.Drawing.Point(243, 37);
            this.txtBoxMaxUV.Name = "txtBoxMaxUV";
            this.txtBoxMaxUV.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMaxUV.TabIndex = 58;
            this.txtBoxMaxUV.Text = "1,5";
            this.txtBoxMaxUV.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // chart
            // 
            chartArea6.Name = "ChartArea1";
            this.chart.ChartAreas.Add(chartArea6);
            legend6.Name = "Legend1";
            this.chart.Legends.Add(legend6);
            this.chart.Location = new System.Drawing.Point(242, 8);
            this.chart.Name = "chart";
            series11.BorderWidth = 3;
            series11.ChartArea = "ChartArea1";
            series11.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series11.Legend = "Legend1";
            series11.MarkerSize = 7;
            series11.Name = "U";
            series11.XValueType = System.Windows.Forms.DataVisualization.Charting.ChartValueType.Double;
            series11.YValueType = System.Windows.Forms.DataVisualization.Charting.ChartValueType.Double;
            series12.BorderWidth = 2;
            series12.ChartArea = "ChartArea1";
            series12.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series12.Legend = "Legend1";
            series12.Name = "U2";
            this.chart.Series.Add(series11);
            this.chart.Series.Add(series12);
            this.chart.Size = new System.Drawing.Size(1031, 650);
            this.chart.TabIndex = 44;
            this.chart.Text = "chart1";
            // 
            // btnSolve
            // 
            this.btnSolve.Location = new System.Drawing.Point(16, 579);
            this.btnSolve.Name = "btnSolve";
            this.btnSolve.Size = new System.Drawing.Size(86, 23);
            this.btnSolve.TabIndex = 45;
            this.btnSolve.Text = "Решить";
            this.btnSolve.UseVisualStyleBackColor = true;
            this.btnSolve.Click += new System.EventHandler(this.btnSolve_Click);
            // 
            // btnPlot
            // 
            this.btnPlot.Enabled = false;
            this.btnPlot.Location = new System.Drawing.Point(16, 608);
            this.btnPlot.Name = "btnPlot";
            this.btnPlot.Size = new System.Drawing.Size(86, 23);
            this.btnPlot.TabIndex = 46;
            this.btnPlot.Text = "Нарисовать";
            this.btnPlot.UseVisualStyleBackColor = true;
            this.btnPlot.Click += new System.EventHandler(this.btnPlot_Click);
            // 
            // trBarT
            // 
            this.trBarT.Location = new System.Drawing.Point(312, 658);
            this.trBarT.Maximum = 49;
            this.trBarT.Name = "trBarT";
            this.trBarT.Size = new System.Drawing.Size(827, 45);
            this.trBarT.TabIndex = 47;
            this.trBarT.Scroll += new System.EventHandler(this.trBarT_Scroll);
            // 
            // prBarSolve
            // 
            this.prBarSolve.Location = new System.Drawing.Point(16, 661);
            this.prBarSolve.Maximum = 4;
            this.prBarSolve.Name = "prBarSolve";
            this.prBarSolve.Size = new System.Drawing.Size(167, 23);
            this.prBarSolve.Step = 1;
            this.prBarSolve.Style = System.Windows.Forms.ProgressBarStyle.Continuous;
            this.prBarSolve.TabIndex = 48;
            // 
            // rdBtnTmr
            // 
            this.rdBtnTmr.AutoSize = true;
            this.rdBtnTmr.Location = new System.Drawing.Point(111, 609);
            this.rdBtnTmr.Name = "rdBtnTmr";
            this.rdBtnTmr.Size = new System.Drawing.Size(86, 17);
            this.rdBtnTmr.TabIndex = 49;
            this.rdBtnTmr.Text = "С таймером";
            this.rdBtnTmr.UseVisualStyleBackColor = true;
            // 
            // btnStopTimer
            // 
            this.btnStopTimer.Location = new System.Drawing.Point(108, 632);
            this.btnStopTimer.Name = "btnStopTimer";
            this.btnStopTimer.Size = new System.Drawing.Size(75, 23);
            this.btnStopTimer.TabIndex = 50;
            this.btnStopTimer.Text = "Stop Timer";
            this.btnStopTimer.UseVisualStyleBackColor = true;
            this.btnStopTimer.Click += new System.EventHandler(this.btnStopTimer_Click);
            // 
            // timerT
            // 
            this.timerT.Interval = 1;
            this.timerT.Tick += new System.EventHandler(this.timerT_Tick);
            // 
            // btnAbout
            // 
            this.btnAbout.Location = new System.Drawing.Point(12, 12);
            this.btnAbout.Name = "btnAbout";
            this.btnAbout.Size = new System.Drawing.Size(75, 23);
            this.btnAbout.TabIndex = 84;
            this.btnAbout.Text = "Об ур-нии";
            this.btnAbout.UseVisualStyleBackColor = true;
            this.btnAbout.Click += new System.EventHandler(this.btnAbout_Click);
            // 
            // propertyGrid1
            // 
            this.propertyGrid1.Location = new System.Drawing.Point(-4, 0);
            this.propertyGrid1.Name = "propertyGrid1";
            this.propertyGrid1.PropertySort = System.Windows.Forms.PropertySort.Alphabetical;
            this.propertyGrid1.Size = new System.Drawing.Size(171, 369);
            this.propertyGrid1.TabIndex = 85;
            this.propertyGrid1.PropertyValueChanged += new System.Windows.Forms.PropertyValueChangedEventHandler(this.propertyGrid_PropertyValueChanged);
            // 
            // tabControl1
            // 
            this.tabControl1.Controls.Add(this.tabPage1);
            this.tabControl1.Controls.Add(this.tabPage2);
            this.tabControl1.Location = new System.Drawing.Point(12, 41);
            this.tabControl1.Name = "tabControl1";
            this.tabControl1.SelectedIndex = 0;
            this.tabControl1.Size = new System.Drawing.Size(171, 395);
            this.tabControl1.TabIndex = 86;
            // 
            // tabPage1
            // 
            this.tabPage1.Controls.Add(this.propertyGrid1);
            this.tabPage1.Location = new System.Drawing.Point(4, 22);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage1.Size = new System.Drawing.Size(163, 369);
            this.tabPage1.TabIndex = 0;
            this.tabPage1.Text = "1-ое ур-ние";
            this.tabPage1.UseVisualStyleBackColor = true;
            // 
            // tabPage2
            // 
            this.tabPage2.Controls.Add(this.propertyGrid2);
            this.tabPage2.Location = new System.Drawing.Point(4, 22);
            this.tabPage2.Name = "tabPage2";
            this.tabPage2.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage2.Size = new System.Drawing.Size(163, 369);
            this.tabPage2.TabIndex = 1;
            this.tabPage2.Text = "2-ое ур-ние";
            this.tabPage2.UseVisualStyleBackColor = true;
            // 
            // propertyGrid2
            // 
            this.propertyGrid2.Location = new System.Drawing.Point(-4, 0);
            this.propertyGrid2.Name = "propertyGrid2";
            this.propertyGrid2.PropertySort = System.Windows.Forms.PropertySort.Alphabetical;
            this.propertyGrid2.Size = new System.Drawing.Size(171, 369);
            this.propertyGrid2.TabIndex = 86;
            this.propertyGrid2.PropertyValueChanged += new System.Windows.Forms.PropertyValueChangedEventHandler(this.propertyGrid_PropertyValueChanged);
            // 
            // WindowPDE
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1284, 711);
            this.Controls.Add(this.tabControl1);
            this.Controls.Add(this.btnAbout);
            this.Controls.Add(this.lblError);
            this.Controls.Add(this.checkBox2ndEq);
            this.Controls.Add(this.lblMinUV);
            this.Controls.Add(this.lblMaxUV);
            this.Controls.Add(this.btnTune);
            this.Controls.Add(this.txtBoxMinUV);
            this.Controls.Add(this.txtBoxMaxUV);
            this.Controls.Add(this.chart);
            this.Controls.Add(this.btnSolve);
            this.Controls.Add(this.btnPlot);
            this.Controls.Add(this.trBarT);
            this.Controls.Add(this.prBarSolve);
            this.Controls.Add(this.rdBtnTmr);
            this.Controls.Add(this.btnStopTimer);
            this.Name = "WindowPDE";
            this.ShowIcon = false;
            this.Text = "WindowPDE";
            this.Load += new System.EventHandler(this.WindowPDE_Load);
            ((System.ComponentModel.ISupportInitialize)(this.chart)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.trBarT)).EndInit();
            this.tabControl1.ResumeLayout(false);
            this.tabPage1.ResumeLayout(false);
            this.tabPage2.ResumeLayout(false);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label lblError;
        private System.Windows.Forms.CheckBox checkBox2ndEq;
        private System.Windows.Forms.Label lblMinUV;
        private System.Windows.Forms.Label lblMaxUV;
        private System.Windows.Forms.Button btnTune;
        private System.Windows.Forms.TextBox txtBoxMinUV;
        private System.Windows.Forms.TextBox txtBoxMaxUV;
        private System.Windows.Forms.DataVisualization.Charting.Chart chart;
        private System.Windows.Forms.Button btnSolve;
        private System.Windows.Forms.Button btnPlot;
        private System.Windows.Forms.TrackBar trBarT;
        private System.Windows.Forms.ProgressBar prBarSolve;
        private System.Windows.Forms.RadioButton rdBtnTmr;
        private System.Windows.Forms.Button btnStopTimer;
        private System.Windows.Forms.Timer timerT;
        private System.Windows.Forms.Button btnAbout;
        private System.Windows.Forms.PropertyGrid propertyGrid1;
        private System.Windows.Forms.TabControl tabControl1;
        private System.Windows.Forms.TabPage tabPage1;
        private System.Windows.Forms.TabPage tabPage2;
        private System.Windows.Forms.PropertyGrid propertyGrid2;
    }
}