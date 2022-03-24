to run all python code, just type: python3.8 xxx.py
all the codes have similar structures, take plot_yaxis.py at /home/jzl/plot_parkinson_paper/field7/ as an example, a part like lines 117-133 reads the data from the data file, then a part like 136-172 plots out the figure
lines 51-68 are used to set up the double/triple x/y axis

To add/modify one figure, add/modify the part that reads in the data, especially the parameter 'readin_file' and the line numbers in 'linecache.getline(readin_file,xxx)', then add/modify the plotting part. If a new double/triple x/y axis ia needed, modify the parameters in the forward and inverse functions between the origianl x/y axis and the new x/y axis, just like lines 51-68 in plot_yaxis.py

As for the names of the figures, a represents VEN1, b represents VEN2, c represents VEN3 and d represents VEN4

Figure 2 in paper 1: 
use plot_field7.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_height_yaxis.eps (and .png), Figure_b_height_yaxis.eps, Figure_c_height_yaxis.eps, Figure_d_height_yaxis.eps

Figure 3 in paper 1:
use plot_field7.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_wind_yaxis.eps (and .png), Figure_b_wind_yaxis.eps, Figure_c_wind_yaxis.eps, Figure_d_wind_yaxis.eps

Figure 4 in paper 1:
use plot_field7.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_temperature_yaxis.eps (and .png), Figure_b_temperature_yaxis.eps, Figure_c_temperature_yaxis.eps, Figure_d_temperature_yaxis.eps

Figure 5 in paper 1:
use new_heat_balance.py at /home/jzl/plot_parkinson_paper/field6/
the output contains Figure_heat_a.eps (and .png), Figure_heat_b.eps, Figure_heat_c.eps, Figure_heat_d.eps


Figure 6a in paper 1:
use plot_new_figures.py at /home/jzl/plot_parkinson_paper/new_field/
the output contains Figure_nir.eps (and .png)

Figure 6b and 6c in paper 1:
use plot_new_figures_20210926.py at /home/jzl/plot_parkinson_paper/new_field_20210921/
the output contains Figure_ven4_0p7_temperature_yaxis.eps (and .png), Figure_ven4_1p3_temperature_yaxis.eps

Figure 2 in paper 2: same as Figure 2 in paper 1

Figure 3 in paper 2: same as Figure 4 in paper 1

Figure 4 in paper 2:
use plot_yaxis.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_o2_yaxis.eps (and .png), Figure_b_o2_yaxis.eps, Figure_c_o2_yaxis.eps, Figure_d_o2_yaxis.eps

Figure 5 in paper 2:
use plot_o2_yaxis.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_o2_plot_yaxis.eps (and .png), Figure_b_o2_plot_yaxis.eps, Figure_c_o2_plot_yaxis.eps, Figure_d_o2_plot_yaxis.eps

Figure 6 in paper 2:
use plot_yaxis.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_o1_yaxis.eps (and .png), Figure_b_o1_yaxis.eps, Figure_c_o1_yaxis.eps, Figure_d_o1_yaxis.eps

Figure 7 in paper 2:
use plot_field5.py at /home/jzl/plot_parkinson_paper/field5/
the output contains Figure_a_o3.eps (and _yaxis.png), Figure_b_o3.eps, Figure_c_o3.eps, Figure_d_o3.eps

Figure 8 in paper 2:
use plot_yaxis.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_co2_yaxis.eps (and .png), Figure_b_co2_yaxis.eps, Figure_c_co2_yaxis.eps, Figure_d_co2_yaxis.eps

Figure 9 in paper 2:
use plot_yaxis.py at /home/jzl/plot_parkinson_paper/field7/
the output contains Figure_a_co_yaxis.eps (and .png), Figure_b_co_yaxis.eps, Figure_c_co_yaxis.eps, Figure_d_co_yaxis.eps

Figure 10 in paper 2:
use plot_field5.py at /home/jzl/plot_parkinson_paper/field5/
the output contains Figure_a_so.eps (and _yaxis.png), Figure_b_so.eps, Figure_c_so.eps, Figure_d_so.eps

Figure 11 in paper 2:
use plot_field5.py at /home/jzl/plot_parkinson_paper/field5/
the output contains Figure_a_so2.eps (and _yaxis.png), Figure_b_so2.eps, Figure_c_so2.eps, Figure_d_so2.eps

Figure 12 in paper 2:
use plot_field6.py at /home/jzl/plot_parkinson_paper/field6/
the output contains Figure_a_noir.eps (and _yaxis.png), Figure_b_noir.eps, Figure_c_noir.eps, Figure_d_noir.eps

Figure 13 in paper 2:
use plot_field5.py at /home/jzl/plot_parkinson_paper/field5/
the output contains Figure_a_o2ir.eps (and _yaxis.png), Figure_b_o2ir.eps, Figure_c_o2ir.eps, Figure_d_o2ir.eps

Figure 14 in paper 2:
use plot_field5.py at /home/jzl/plot_parkinson_paper/field5/
the output contains Figure_a_ohir.eps (and _yaxis.png), Figure_b_ohir.eps, Figure_c_ohir.eps, Figure_d_ohir.eps




