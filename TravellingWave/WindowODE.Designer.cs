namespace TravellingWave
{
    partial class WindowODE
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
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend1 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series2 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series3 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series4 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend2 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series5 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series6 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series7 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series8 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series9 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series10 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.txtBoxMaxUVT = new System.Windows.Forms.TextBox();
            this.rdBtnPlotAll = new System.Windows.Forms.RadioButton();
            this.txtBoxMinUVT = new System.Windows.Forms.TextBox();
            this.btnTuneT = new System.Windows.Forms.Button();
            this.rdBtnTmr = new System.Windows.Forms.RadioButton();
            this.lblMaxPhase = new System.Windows.Forms.Label();
            this.txtBoxMaxUVPhase = new System.Windows.Forms.TextBox();
            this.lblMinPhase = new System.Windows.Forms.Label();
            this.txtBoxMinUVPhase = new System.Windows.Forms.TextBox();
            this.btnTunePhase = new System.Windows.Forms.Button();
            this.lblError = new System.Windows.Forms.Label();
            this.chart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.trBarT = new System.Windows.Forms.TrackBar();
            this.chartPhase = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.prBarSolve = new System.Windows.Forms.ProgressBar();
            this.lblMinUV = new System.Windows.Forms.Label();
            this.panel1 = new System.Windows.Forms.Panel();
            this.btnStopTimer = new System.Windows.Forms.Button();
            this.lblMaxUV = new System.Windows.Forms.Label();
            this.btnPlot = new System.Windows.Forms.Button();
            this.btnSolve = new System.Windows.Forms.Button();
            this.timerT = new System.Windows.Forms.Timer(this.components);
            this.tabControl1 = new System.Windows.Forms.TabControl();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.propertyGrid1 = new System.Windows.Forms.PropertyGrid();
            this.tabPage2 = new System.Windows.Forms.TabPage();
            this.propertyGrid2 = new System.Windows.Forms.PropertyGrid();
            this.checkBox2ndEq = new System.Windows.Forms.CheckBox();
            this.btnAbout = new System.Windows.Forms.Button();
            this.btnTickFaster = new System.Windows.Forms.Button();
            this.btnTickSlower = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.chart)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.trBarT)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.chartPhase)).BeginInit();
            this.panel1.SuspendLayout();
            this.tabControl1.SuspendLayout();
            this.tabPage1.SuspendLayout();
            this.tabPage2.SuspendLayout();
            this.SuspendLayout();
            // 
            // txtBoxMaxUVT
            // 
            this.txtBoxMaxUVT.Location = new System.Drawing.Point(247, 25);
            this.txtBoxMaxUVT.Name = "txtBoxMaxUVT";
            this.txtBoxMaxUVT.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMaxUVT.TabIndex = 79;
            this.txtBoxMaxUVT.Text = "2,0";
            this.txtBoxMaxUVT.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // rdBtnPlotAll
            // 
            this.rdBtnPlotAll.AutoSize = true;
            this.rdBtnPlotAll.Checked = true;
            this.rdBtnPlotAll.Location = new System.Drawing.Point(3, 26);
            this.rdBtnPlotAll.Name = "rdBtnPlotAll";
            this.rdBtnPlotAll.Size = new System.Drawing.Size(57, 17);
            this.rdBtnPlotAll.TabIndex = 1;
            this.rdBtnPlotAll.TabStop = true;
            this.rdBtnPlotAll.Text = "Plot All";
            this.rdBtnPlotAll.UseVisualStyleBackColor = true;
            // 
            // txtBoxMinUVT
            // 
            this.txtBoxMinUVT.Location = new System.Drawing.Point(247, 276);
            this.txtBoxMinUVT.Name = "txtBoxMinUVT";
            this.txtBoxMinUVT.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMinUVT.TabIndex = 82;
            this.txtBoxMinUVT.Text = "-2,1";
            this.txtBoxMinUVT.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // btnTuneT
            // 
            this.btnTuneT.Location = new System.Drawing.Point(247, 152);
            this.btnTuneT.Name = "btnTuneT";
            this.btnTuneT.Size = new System.Drawing.Size(42, 23);
            this.btnTuneT.TabIndex = 71;
            this.btnTuneT.Text = "Zoom";
            this.btnTuneT.UseVisualStyleBackColor = true;
            this.btnTuneT.Click += new System.EventHandler(this.btnTuneT_Click);
            // 
            // rdBtnTmr
            // 
            this.rdBtnTmr.AutoSize = true;
            this.rdBtnTmr.Location = new System.Drawing.Point(3, 3);
            this.rdBtnTmr.Name = "rdBtnTmr";
            this.rdBtnTmr.Size = new System.Drawing.Size(76, 17);
            this.rdBtnTmr.TabIndex = 0;
            this.rdBtnTmr.Text = "With Timer";
            this.rdBtnTmr.UseVisualStyleBackColor = true;
            // 
            // lblMaxPhase
            // 
            this.lblMaxPhase.AutoSize = true;
            this.lblMaxPhase.Location = new System.Drawing.Point(203, 365);
            this.lblMaxPhase.Name = "lblMaxPhase";
            this.lblMaxPhase.Size = new System.Drawing.Size(44, 13);
            this.lblMaxPhase.TabIndex = 83;
            this.lblMaxPhase.Text = "max u,v";
            // 
            // txtBoxMaxUVPhase
            // 
            this.txtBoxMaxUVPhase.Location = new System.Drawing.Point(247, 362);
            this.txtBoxMaxUVPhase.Name = "txtBoxMaxUVPhase";
            this.txtBoxMaxUVPhase.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMaxUVPhase.TabIndex = 77;
            this.txtBoxMaxUVPhase.Text = "2,0";
            this.txtBoxMaxUVPhase.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // lblMinPhase
            // 
            this.lblMinPhase.AutoSize = true;
            this.lblMinPhase.Location = new System.Drawing.Point(203, 602);
            this.lblMinPhase.Name = "lblMinPhase";
            this.lblMinPhase.Size = new System.Drawing.Size(41, 13);
            this.lblMinPhase.TabIndex = 84;
            this.lblMinPhase.Text = "min u,v";
            // 
            // txtBoxMinUVPhase
            // 
            this.txtBoxMinUVPhase.Location = new System.Drawing.Point(247, 597);
            this.txtBoxMinUVPhase.Name = "txtBoxMinUVPhase";
            this.txtBoxMinUVPhase.Size = new System.Drawing.Size(34, 20);
            this.txtBoxMinUVPhase.TabIndex = 81;
            this.txtBoxMinUVPhase.Text = "-1,0";
            this.txtBoxMinUVPhase.TextAlign = System.Windows.Forms.HorizontalAlignment.Center;
            // 
            // btnTunePhase
            // 
            this.btnTunePhase.Location = new System.Drawing.Point(247, 489);
            this.btnTunePhase.Name = "btnTunePhase";
            this.btnTunePhase.Size = new System.Drawing.Size(42, 23);
            this.btnTunePhase.TabIndex = 75;
            this.btnTunePhase.Text = "Zoom";
            this.btnTunePhase.UseVisualStyleBackColor = true;
            this.btnTunePhase.Click += new System.EventHandler(this.btnTunePhase_Click);
            // 
            // lblError
            // 
            this.lblError.AutoSize = true;
            this.lblError.Font = new System.Drawing.Font("Microsoft Sans Serif", 9.75F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.lblError.ForeColor = System.Drawing.Color.DarkRed;
            this.lblError.Location = new System.Drawing.Point(193, 665);
            this.lblError.Name = "lblError";
            this.lblError.Size = new System.Drawing.Size(102, 16);
            this.lblError.TabIndex = 85;
            this.lblError.Text = "Error occured";
            this.lblError.Visible = false;
            // 
            // chart
            // 
            chartArea1.Name = "ChartArea1";
            this.chart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.chart.Legends.Add(legend1);
            this.chart.Location = new System.Drawing.Point(247, 8);
            this.chart.Name = "chart";
            series1.BorderWidth = 3;
            series1.ChartArea = "ChartArea1";
            series1.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series1.Legend = "Legend1";
            series1.Name = "U";
            series2.BorderWidth = 3;
            series2.ChartArea = "ChartArea1";
            series2.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series2.Legend = "Legend1";
            series2.Name = "V";
            series3.BorderWidth = 2;
            series3.ChartArea = "ChartArea1";
            series3.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series3.Legend = "Legend1";
            series3.Name = "U2";
            series4.BorderWidth = 2;
            series4.ChartArea = "ChartArea1";
            series4.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series4.Legend = "Legend1";
            series4.Name = "V2";
            this.chart.Series.Add(series1);
            this.chart.Series.Add(series2);
            this.chart.Series.Add(series3);
            this.chart.Series.Add(series4);
            this.chart.Size = new System.Drawing.Size(1030, 325);
            this.chart.TabIndex = 72;
            this.chart.Text = "chart1";
            // 
            // trBarT
            // 
            this.trBarT.Location = new System.Drawing.Point(328, 658);
            this.trBarT.Name = "trBarT";
            this.trBarT.Size = new System.Drawing.Size(788, 45);
            this.trBarT.TabIndex = 70;
            this.trBarT.Scroll += new System.EventHandler(this.trBarT_Scroll);
            // 
            // chartPhase
            // 
            chartArea2.Name = "ChartArea1";
            this.chartPhase.ChartAreas.Add(chartArea2);
            legend2.Name = "Legend1";
            this.chartPhase.Legends.Add(legend2);
            this.chartPhase.Location = new System.Drawing.Point(247, 333);
            this.chartPhase.Name = "chartPhase";
            series5.BorderWidth = 5;
            series5.ChartArea = "ChartArea1";
            series5.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series5.Legend = "Legend1";
            series5.Name = "V(U)";
            series6.BorderWidth = 3;
            series6.ChartArea = "ChartArea1";
            series6.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series6.Legend = "Legend1";
            series6.Name = "1 nс V(U)";
            series7.BorderWidth = 3;
            series7.ChartArea = "ChartArea1";
            series7.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series7.Legend = "Legend1";
            series7.Name = "2 nc V(U)";
            series8.BorderWidth = 2;
            series8.ChartArea = "ChartArea1";
            series8.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series8.Legend = "Legend1";
            series8.Name = "V2";
            series8.Palette = System.Windows.Forms.DataVisualization.Charting.ChartColorPalette.BrightPastel;
            series9.BorderWidth = 2;
            series9.ChartArea = "ChartArea1";
            series9.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series9.Legend = "Legend1";
            series9.Name = "1 nc V2";
            series9.Palette = System.Windows.Forms.DataVisualization.Charting.ChartColorPalette.BrightPastel;
            series10.BorderWidth = 2;
            series10.ChartArea = "ChartArea1";
            series10.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            series10.Legend = "Legend1";
            series10.Name = "2 nc V2";
            series10.Palette = System.Windows.Forms.DataVisualization.Charting.ChartColorPalette.BrightPastel;
            this.chartPhase.Series.Add(series5);
            this.chartPhase.Series.Add(series6);
            this.chartPhase.Series.Add(series7);
            this.chartPhase.Series.Add(series8);
            this.chartPhase.Series.Add(series9);
            this.chartPhase.Series.Add(series10);
            this.chartPhase.Size = new System.Drawing.Size(1030, 325);
            this.chartPhase.TabIndex = 74;
            this.chartPhase.Text = "chart2";
            // 
            // prBarSolve
            // 
            this.prBarSolve.Location = new System.Drawing.Point(31, 661);
            this.prBarSolve.Maximum = 3;
            this.prBarSolve.Name = "prBarSolve";
            this.prBarSolve.Size = new System.Drawing.Size(156, 23);
            this.prBarSolve.TabIndex = 47;
            // 
            // lblMinUV
            // 
            this.lblMinUV.AutoSize = true;
            this.lblMinUV.Location = new System.Drawing.Point(206, 279);
            this.lblMinUV.Name = "lblMinUV";
            this.lblMinUV.Size = new System.Drawing.Size(41, 13);
            this.lblMinUV.TabIndex = 76;
            this.lblMinUV.Text = "min u,v";
            // 
            // panel1
            // 
            this.panel1.Controls.Add(this.rdBtnPlotAll);
            this.panel1.Controls.Add(this.rdBtnTmr);
            this.panel1.Location = new System.Drawing.Point(115, 579);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(80, 47);
            this.panel1.TabIndex = 69;
            // 
            // btnStopTimer
            // 
            this.btnStopTimer.Location = new System.Drawing.Point(112, 632);
            this.btnStopTimer.Name = "btnStopTimer";
            this.btnStopTimer.Size = new System.Drawing.Size(75, 23);
            this.btnStopTimer.TabIndex = 68;
            this.btnStopTimer.Text = "Stop Timer";
            this.btnStopTimer.UseVisualStyleBackColor = true;
            this.btnStopTimer.Click += new System.EventHandler(this.btnStopTimer_Click);
            // 
            // lblMaxUV
            // 
            this.lblMaxUV.AutoSize = true;
            this.lblMaxUV.Location = new System.Drawing.Point(203, 28);
            this.lblMaxUV.Name = "lblMaxUV";
            this.lblMaxUV.Size = new System.Drawing.Size(44, 13);
            this.lblMaxUV.TabIndex = 73;
            this.lblMaxUV.Text = "max u,v";
            // 
            // btnPlot
            // 
            this.btnPlot.Enabled = false;
            this.btnPlot.Location = new System.Drawing.Point(31, 608);
            this.btnPlot.Name = "btnPlot";
            this.btnPlot.Size = new System.Drawing.Size(75, 23);
            this.btnPlot.TabIndex = 67;
            this.btnPlot.Text = "Plot";
            this.btnPlot.UseVisualStyleBackColor = true;
            this.btnPlot.Click += new System.EventHandler(this.btnPlot_Click);
            // 
            // btnSolve
            // 
            this.btnSolve.Location = new System.Drawing.Point(31, 579);
            this.btnSolve.Name = "btnSolve";
            this.btnSolve.Size = new System.Drawing.Size(75, 23);
            this.btnSolve.TabIndex = 66;
            this.btnSolve.Text = "Solve";
            this.btnSolve.UseVisualStyleBackColor = true;
            this.btnSolve.Click += new System.EventHandler(this.btnSolve_Click);
            // 
            // timerT
            // 
            this.timerT.Interval = 50;
            this.timerT.Tick += new System.EventHandler(this.timerT_Tick);
            // 
            // tabControl1
            // 
            this.tabControl1.Controls.Add(this.tabPage1);
            this.tabControl1.Controls.Add(this.tabPage2);
            this.tabControl1.Location = new System.Drawing.Point(8, 43);
            this.tabControl1.Name = "tabControl1";
            this.tabControl1.SelectedIndex = 0;
            this.tabControl1.Size = new System.Drawing.Size(176, 328);
            this.tabControl1.TabIndex = 86;
            // 
            // tabPage1
            // 
            this.tabPage1.Controls.Add(this.propertyGrid1);
            this.tabPage1.Location = new System.Drawing.Point(4, 22);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage1.Size = new System.Drawing.Size(168, 302);
            this.tabPage1.TabIndex = 0;
            this.tabPage1.Text = "First eq";
            this.tabPage1.UseVisualStyleBackColor = true;
            // 
            // propertyGrid1
            // 
            this.propertyGrid1.Location = new System.Drawing.Point(0, 0);
            this.propertyGrid1.Name = "propertyGrid1";
            this.propertyGrid1.PropertySort = System.Windows.Forms.PropertySort.Alphabetical;
            this.propertyGrid1.Size = new System.Drawing.Size(168, 300);
            this.propertyGrid1.TabIndex = 0;
            this.propertyGrid1.SelectedGridItemChanged += new System.Windows.Forms.SelectedGridItemChangedEventHandler(this.propertyGrid1_SelectedGridItemChanged);
            // 
            // tabPage2
            // 
            this.tabPage2.Controls.Add(this.propertyGrid2);
            this.tabPage2.Location = new System.Drawing.Point(4, 22);
            this.tabPage2.Name = "tabPage2";
            this.tabPage2.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage2.Size = new System.Drawing.Size(168, 302);
            this.tabPage2.TabIndex = 1;
            this.tabPage2.Text = "Second eq";
            this.tabPage2.UseVisualStyleBackColor = true;
            // 
            // propertyGrid2
            // 
            this.propertyGrid2.Location = new System.Drawing.Point(0, 0);
            this.propertyGrid2.Name = "propertyGrid2";
            this.propertyGrid2.PropertySort = System.Windows.Forms.PropertySort.Alphabetical;
            this.propertyGrid2.Size = new System.Drawing.Size(168, 300);
            this.propertyGrid2.TabIndex = 1;
            this.propertyGrid2.SelectedGridItemChanged += new System.Windows.Forms.SelectedGridItemChangedEventHandler(this.propertyGrid2_SelectedGridItemChanged);
            // 
            // checkBox2ndEq
            // 
            this.checkBox2ndEq.AutoSize = true;
            this.checkBox2ndEq.Location = new System.Drawing.Point(12, 377);
            this.checkBox2ndEq.Name = "checkBox2ndEq";
            this.checkBox2ndEq.Size = new System.Drawing.Size(59, 17);
            this.checkBox2ndEq.TabIndex = 1;
            this.checkBox2ndEq.Text = "2nd eq";
            this.checkBox2ndEq.UseVisualStyleBackColor = true;
            this.checkBox2ndEq.CheckedChanged += new System.EventHandler(this.checkBox2ndEq_CheckedChanged);
            // 
            // btnAbout
            // 
            this.btnAbout.Location = new System.Drawing.Point(8, 12);
            this.btnAbout.Name = "btnAbout";
            this.btnAbout.Size = new System.Drawing.Size(75, 23);
            this.btnAbout.TabIndex = 87;
            this.btnAbout.Text = "About";
            this.btnAbout.UseVisualStyleBackColor = true;
            this.btnAbout.Click += new System.EventHandler(this.btnAbout_Click);
            // 
            // btnTickFaster
            // 
            this.btnTickFaster.Location = new System.Drawing.Point(147, 553);
            this.btnTickFaster.Name = "btnTickFaster";
            this.btnTickFaster.Size = new System.Drawing.Size(26, 23);
            this.btnTickFaster.TabIndex = 90;
            this.btnTickFaster.Text = "+";
            this.btnTickFaster.UseVisualStyleBackColor = true;
            this.btnTickFaster.Click += new System.EventHandler(this.btnTickFaster_Click);
            // 
            // btnTickSlower
            // 
            this.btnTickSlower.Location = new System.Drawing.Point(115, 553);
            this.btnTickSlower.Name = "btnTickSlower";
            this.btnTickSlower.Size = new System.Drawing.Size(26, 23);
            this.btnTickSlower.TabIndex = 89;
            this.btnTickSlower.Text = "-";
            this.btnTickSlower.UseVisualStyleBackColor = true;
            this.btnTickSlower.Click += new System.EventHandler(this.btnTickSlower_Click);
            // 
            // WindowODE
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1284, 711);
            this.Controls.Add(this.btnTickFaster);
            this.Controls.Add(this.btnTickSlower);
            this.Controls.Add(this.btnAbout);
            this.Controls.Add(this.checkBox2ndEq);
            this.Controls.Add(this.tabControl1);
            this.Controls.Add(this.txtBoxMaxUVT);
            this.Controls.Add(this.txtBoxMinUVT);
            this.Controls.Add(this.btnTuneT);
            this.Controls.Add(this.lblMaxPhase);
            this.Controls.Add(this.txtBoxMaxUVPhase);
            this.Controls.Add(this.lblMinPhase);
            this.Controls.Add(this.txtBoxMinUVPhase);
            this.Controls.Add(this.btnTunePhase);
            this.Controls.Add(this.lblError);
            this.Controls.Add(this.chart);
            this.Controls.Add(this.trBarT);
            this.Controls.Add(this.chartPhase);
            this.Controls.Add(this.prBarSolve);
            this.Controls.Add(this.lblMinUV);
            this.Controls.Add(this.panel1);
            this.Controls.Add(this.btnStopTimer);
            this.Controls.Add(this.lblMaxUV);
            this.Controls.Add(this.btnPlot);
            this.Controls.Add(this.btnSolve);
            this.Name = "WindowODE";
            this.ShowIcon = false;
            this.Text = "WindowODE";
            this.Load += new System.EventHandler(this.WindowODE_Load);
            ((System.ComponentModel.ISupportInitialize)(this.chart)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.trBarT)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.chartPhase)).EndInit();
            this.panel1.ResumeLayout(false);
            this.panel1.PerformLayout();
            this.tabControl1.ResumeLayout(false);
            this.tabPage1.ResumeLayout(false);
            this.tabPage2.ResumeLayout(false);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.TextBox txtBoxMaxUVT;
        private System.Windows.Forms.RadioButton rdBtnPlotAll;
        private System.Windows.Forms.TextBox txtBoxMinUVT;
        private System.Windows.Forms.Button btnTuneT;
        private System.Windows.Forms.RadioButton rdBtnTmr;
        private System.Windows.Forms.Label lblMaxPhase;
        private System.Windows.Forms.TextBox txtBoxMaxUVPhase;
        private System.Windows.Forms.Label lblMinPhase;
        private System.Windows.Forms.TextBox txtBoxMinUVPhase;
        private System.Windows.Forms.Button btnTunePhase;
        private System.Windows.Forms.Label lblError;
        private System.Windows.Forms.DataVisualization.Charting.Chart chart;
        private System.Windows.Forms.TrackBar trBarT;
        private System.Windows.Forms.DataVisualization.Charting.Chart chartPhase;
        private System.Windows.Forms.ProgressBar prBarSolve;
        private System.Windows.Forms.Label lblMinUV;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Button btnStopTimer;
        private System.Windows.Forms.Label lblMaxUV;
        private System.Windows.Forms.Button btnPlot;
        private System.Windows.Forms.Button btnSolve;
        private System.Windows.Forms.Timer timerT;
        private System.Windows.Forms.TabControl tabControl1;
        private System.Windows.Forms.TabPage tabPage1;
        private System.Windows.Forms.PropertyGrid propertyGrid1;
        private System.Windows.Forms.TabPage tabPage2;
        private System.Windows.Forms.PropertyGrid propertyGrid2;
        private System.Windows.Forms.CheckBox checkBox2ndEq;
        private System.Windows.Forms.Button btnAbout;
        private System.Windows.Forms.Button btnTickFaster;
        private System.Windows.Forms.Button btnTickSlower;
    }
}