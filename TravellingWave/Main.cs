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
    public partial class Main : Form
    {
        public Main()
        {
            InitializeComponent();
        }

        private void btnLoadODE_Click(object sender, EventArgs e)
        {
            WindowODE o = new WindowODE();
            o.ShowDialog();
        }

        private void btnLoadPDE_Click(object sender, EventArgs e)
        {
            WindowPDE p = new WindowPDE();
            p.ShowDialog();
        }
    }
}
