//Input parameters..............................
//display: 0:Fit + plots are done,  1:only plots are done using the last found parameters
//m_0: Lower limit of the fit range
//m_cut: Upper limit of the fit range
//t1: t < t1 Early type galaxies (E/S0)
//t2: t1 < t < t2 Spirals
//n_gal_fit: Number of galaxies used for the fit (write .0). They will be printed in mzt_small.txt

#define display 0
#define m_0 19.0
#define m_cut 24.0
#define t1 16.0
#define t2 54.0
#define n_gal_fit 10000.0

#define bins_m 60 //Magnitude bins in the prior plot
#define bins_z 60 //Redshift bins in the prior plot
#define bins_t 3 //Type bins in the prior plot

//File Paths......................................
//INPUT:
//cat_file: a list of values: observed mag in certain filt, z_true, template index for all the
//galaxies of the mock catalog. This list provides the necessary data to fit the prior function.
char cat_file[1000] = "../../Data/Catalog/mock.r260.n1e6.s10.121027_mzt.txt";

//OUTPUT:
//prior_para_file: file with the resulting parameters of the fit
//prior_cat_file: If we don't want to use all the point in mzt.txt to do the calibration we can select a fraction of
//them randomly. They will be printed in this file mzt_small.txt.
//prior_py_file: A file that also contains the values of the parameters but in a python format to directly copy&paste 
//in the BPZ prior_{}.py script
char prior_para_file[1000] = "../../Data/Photoz/PAU/PRIOR/prior_para.txt";
char small_cat_file[1000] = "../../Data/Photoz/PAU/PRIOR/mzt_small.txt";
char prior_py_file[1000] = "../../Data/Photoz/PAU/PRIOR/py_prior.txt";
char plot_file[1000] = "../../Data/Photoz/PAU/PRIOR/plots/";
char plot_txt_file[1000] = "../../Data/Photoz/PAU/PRIOR/txt/";
