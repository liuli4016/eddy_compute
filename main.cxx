#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <math.h>
#include <iterator>
#include <cmath>
#include <iostream>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <ctime>
#include <netcdf.h>

#include <vtkVersion.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPointSource.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkExtractSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelectedBlock.h>
#include <vtkMarchingCubes.h>
#include <vtkInformation.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataWriter.h>
#include <vtkImplicitModeller.h>
#include <vtkCell.h>
#include "vtkStructuredGrid.h"
#include <vtkXMLStructuredGridWriter.h>
#include "vtkDebugLeaks.h"
#include <vtkDiscreteMarchingCubes.h>
#include "vtkStructuredPointsReader.h"
#include "vtkPolyDataConnectivityFilter.h"
#include <vtkFillHolesFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkReverseSense.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include <vtkDataArray.h>
#include <vtkContourFilter.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkType.h>

#include "FeatureTrack.h"
#include "input.h"
#include "classifier.h"


using namespace std;

int const Consts::MAXSPLIT = 500;
float const Consts::RADIUS   = 1;
float const Consts::FRAMEDIST = 10.0;
float const Consts::NODEDIST = 2.0;
long  const Consts::MAXOBJS = 100000;
int  const Consts::MAXLEN = 10000;
int const Consts::TRAKTABLE_MAXLEN = 409600;
float const Consts::DEFAULT_TOLERANCE = 0.0;
bool const Consts::YES = true;
bool const Consts::NO = false;
int const Consts::DEFAULT_TIMESTEP_NUM = 1200;
bool const Consts::WriteTrackInfo = true;


#define PI 3.14159265

#define DATASET_NAME_x "/grid/x"
#define DATASET_NAME_y "/grid/y"
#define DATASET_NAME_z "/grid/z"
#define DATASET_NAME_swirl "/postprocessing/swirl"
#define DATASET_NAME_T "/flow/T"
#define DATASET_NAME_rho "/flow/rho"
#define DATASET_NAME_u "/flow/u"
#define DATASET_NAME_v "/flow/v"
#define DATASET_NAME_w "/flow/w"









string int2str(int n)
{
   stringstream ss;        
   ss << n; 
   return ss.str();     
}




string getmonth(int i)
{
   string month = "";
	  
   if (i<=31)
      month = "01";
        	 
   else if (i>31 && i<=59)
      month = "02";
      
   else if (i>59 && i<=90)
      month = "03";
      
   else if (i>90 && i<=120)
      month = "04";
      
   else if (i>120 && i<=151)
      month = "05";
      
   else if (i>151 && i<=181)
      month = "06";
      
   else if (i>181 && i<=212)
      month = "07";
      
   else if (i>212 && i<=243)
      month = "08";
      
   else if (i>243 && i<=273)
      month = "09";
      
   else if (i>273 && i<=304)
      month = "10";
      
   else if (i>304 && i<=334)
      month = "11";
      
   else if (i>334 && i<=365)
      month = "12";   
  
   return month;      	 	 
}


string getday(int i)
{
	 string day = "";
	 
	 if (i<=9)
      day = "0" + int2str(i);        	 	  
   else
      day = int2str(i);
      
   return day;
}


string computedate(int i)
{
	 int dumi = 0, dumi2 = 0;
	 string date = "", year = "", month = "", day = "";
        
        if (i<=365)
        {
        	 year = "2005";
        	 dumi = i;
        }
        
        else if (i>365 && i<=730)
        {
        	 year = "2006";
        	 dumi = i-365;
        }
        
        else if (i>730 && i<=1000)
        {
        	 year = "2007";
        	 dumi = i-730;
        }
        	 
        month = getmonth(dumi);
        

        
        if (month == "01")
        	 dumi2 = dumi;
        	 
        else if (month == "02")
        	 dumi2 = dumi-31; 
        	 
        else if (month == "03")
        	 dumi2 = dumi-59;	 	 
        	 
        else if (month == "04")
        	 dumi2 = dumi-90;	
        	 
        else if (month == "05")
        	 dumi2 = dumi-120; 
        
        else if (month == "06")
        	 dumi2 = dumi-151; 	 
        	 
        else if (month == "07")
        	 dumi2 = dumi-181;
        	 
        else if (month == "08")
        	 dumi2 = dumi-212;
        	 
        else if (month == "09")
        	 dumi2 = dumi-243;
        	 
        else if (month == "10")
        	 dumi2 = dumi-273;
        	 
        else if (month == "11")
        	 dumi2 = dumi-304;
        	 
        else if (month == "12")
        	 dumi2 = dumi-334;
        	 
        day = getday(dumi2);
        
        date = year + "/" + month + "/" + day;
        
        return date;
}




float search_zeta(string zetafile, float cx, float cy, float lx, float ly, float ux, float uy)
{
    float cx1, cy1, zeta;
    ifstream fp;
    fp.open(zetafile.c_str(), ios_base::in);
    	
    float radius = 0.5*sqrt((ux-lx)*(ux-lx) + (uy-ly)*(uy-ly));
    float total_zeta = 0;
    long count = 0;

    while (!fp.eof())
    {
        fp >> cx1 >> cy1 >> zeta;

        float dist = sqrt((cx1-cx)*(cx1-cx) + (cy1-cy)*(cy1-cy));

        if ( dist < radius )
        {
           total_zeta = total_zeta + zeta;
           count++;
        }
    }       

    fp.close();

    return total_zeta/count;
}







string precision_time(int const time,int const precision)
{
    int quotient = time, remainder=0;
    int count = 1;
    // THE TOTAL NO OF DIGITS SHOULD EQUAL ATLEAST PRECISION
    while(quotient!=0)
    {
        quotient = quotient/10;
        if(quotient!=0)
        {
            count = count+1;
        }
    }
    string buf1;
    
    for(int i = 0; i< (precision - count);i++)
    {
        buf1+="0";
        
    }
    buf1+=TtoS<int>(time);
    return buf1;
}



vtkDataSet * ReadNcDataFile(string output_base, string FileName, int x_dim, int y_dim, int z_dim, float *threshold1)
{   
    //This function first reads a specific netCDF file and creates a (curvilinear) vtkDataSet from the read data, then returns the content of the file in the vtkDataSet format.

       

    size_t dotplace = FileName.find_last_of(".\\");
    size_t slashplace = FileName.find_last_of("/\\");    
    string GRID_FILE_NAME = FileName.substr(0, slashplace+1) + "NWA_grd.nc";

    string tempstr1 = FileName.substr(slashplace+1, FileName.length());
    string tempstr2 = tempstr1.substr(0, tempstr1.find_last_of(".\\"));



    string zetafile = output_base + tempstr2 + ".zeta";
    
    FILE *fpout1;
    
    if ((fpout1 = fopen(zetafile.c_str(), "w"))==NULL)
       cout << "cannot open zetafile file to write\n";

    
    static int NDIMS = 4;
    int xi_rho = 722;
    int eta_rho = 362;
    //int s_w = 41;
    int s_r = 40;
    int ocean_time = 1;
    double Tcline=0.1;
    double * tempVals = new double [s_r*eta_rho*xi_rho];
    double * saltVals = new double [s_r*eta_rho*xi_rho];
    int retval;
    double * Cs_r = new double [s_r];
   // double * Cs_w = new double [s_w];  
    
    double * pm = new double[eta_rho*xi_rho];
    double * pn = new double[eta_rho*xi_rho];
    double * h = new double [eta_rho*xi_rho];
    double * lon_rho = new double [eta_rho*xi_rho];
    double * lat_rho = new double [eta_rho*xi_rho];
    //double * s_wVals   = new double [s_w];
    double * s_rVals   = new double [s_r];
    double * ZetaVals   = new double [eta_rho*xi_rho];
    double *u_Vals = new double[(xi_rho-1)*eta_rho*s_r];
    double *v_Vals = new double[xi_rho*(eta_rho-1)*s_r];
    //double * ocean_time_Vals = new double [1];
    
    int ncid, pres_varid, temp_varid, salt_varid;
    int lat_varid, lon_varid, Tcline_varid, h_varid, Cs_r_varid, Cs_w_varid, u_Vals_varid, v_Vals_varid, ocean_time__varid, pm_varid, pn_varid ;
    size_t start[NDIMS], count[NDIMS];
    double hc;

    hc= Tcline;

    count[0] = 1;
    count[1] = s_r;
    count[2] = eta_rho;
    count[3] = xi_rho;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    start[0] = 0;
    static size_t start1[] = {0, 0, 0}; /* start at first value */
    static size_t count1[] = {1, eta_rho, xi_rho};


    
    ///-------read grid file first ----------
    if ((retval = nc_open(GRID_FILE_NAME.c_str(), NC_NOWRITE, &ncid)))
        cout<<"couldnt open the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_inq_varid(ncid, "lon_rho", &lat_varid)))
        cout<<"couldnt open 1 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    if ((retval = nc_inq_varid(ncid, "lat_rho", &lon_varid)))
        cout<<"couldnt open 2 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    if ((retval = nc_inq_varid(ncid, "h", &h_varid)))
        cout<<"couldnt open 2.1 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    if ((retval = nc_inq_varid(ncid, "Tcline", &Tcline_varid)))
        cout<<"couldnt open 2.2 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    if ((retval = nc_inq_varid(ncid, "Cs_r", &Cs_r_varid)))
        cout<<"couldnt open 2.3 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_inq_varid(ncid, "pm", &pm_varid)))
        cout<<"couldnt open 2.3.1 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_inq_varid(ncid, "pn", &pn_varid)))
        cout<<"couldnt open 2.3.1 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    
    //if ((retval = nc_inq_varid(ncid, "Cs_w", &Cs_w_varid)))
    //    cout<<"couldnt open 2.4 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    /* Read the coordinate variable data. */

    if ((retval = nc_get_var_double(ncid, pm_varid, pm)))
        cout<<"couldnt open 3 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;

    
    if ((retval = nc_get_var_double(ncid, pn_varid, pn)))
        cout<<"couldnt open 3 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    
    
    
    if ((retval = nc_get_var_double(ncid, lat_varid, lon_rho)))
        cout<<"couldnt open 3 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, lon_varid, lat_rho)))
        cout<<"couldnt open 4 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, h_varid, h)))
        cout<<"couldnt open 4.1 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, Tcline_varid, &Tcline)))
        cout<<"couldnt open 4.2 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, Cs_r_varid, Cs_r)))
        cout<<"couldnt open 4.3 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
   // if ((retval = nc_get_var_double(ncid, Cs_w_varid, Cs_w)))
   //     cout<<"couldnt open 4.4 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;
    
    
    int status = nc_close(ncid);       /* close netCDF dataset */
    if (status != NC_NOERR) cout<<"couldnt open 5 the file :[ "<<GRID_FILE_NAME.c_str()<<" ---- !!!"<<endl;

     

    
    ///-------read the data file now ----------
    
    if ((retval = nc_open(FileName.c_str(), NC_NOWRITE, &ncid)))
        cout<<"couldnt open 6 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
  //  if ((retval = nc_inq_varid(ncid, "s_w", &lon_varid)))
  //      cout<<"couldnt open 7 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
 //   if ((retval = nc_get_var_double(ncid, lon_varid, s_wVals)))
 //       cout<<"couldnt open 8 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;

    if ((retval = nc_inq_varid(ncid, "s_rho", &lon_varid)))
        cout<<"couldnt open 8.1 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, lon_varid, s_rVals)))
        cout<<"couldnt open 8.2 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;

    if ((retval = nc_inq_varid(ncid, "u", & u_Vals_varid)))
        cout<<"couldnt open 8.3 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, u_Vals_varid, u_Vals)))
        cout<<"couldnt open 8.4 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_inq_varid(ncid, "v", & v_Vals_varid)))
        cout<<"couldnt open 8.5 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_var_double(ncid, v_Vals_varid, v_Vals)))
        cout<<"couldnt open 8.6 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    

    if ((retval = nc_inq_varid(ncid, "temp", &temp_varid)))
        cout<<"couldnt open 9.0 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_vara_double(ncid, temp_varid, start, count, tempVals)))
        cout<<"couldnt open 9.1 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;


    if ((retval = nc_inq_varid(ncid, "salt", &salt_varid)))
        cout<<"couldnt open 9.2 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_vara_double(ncid, salt_varid, start, count, saltVals)))
        cout<<"couldnt open 9.3 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;

  
    if ((retval = nc_inq_varid(ncid, "zeta", &lon_varid)))
        cout<<"couldnt open 9.5 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    if ((retval = nc_get_vara_double(ncid, lon_varid, start1, count1, ZetaVals)))
        cout<<"couldnt open 9.7 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;

    
    status = nc_close(ncid);
    if (status != NC_NOERR) cout<<"couldnt open 10 the file :[ "<<FileName.c_str()<<" ---- !!!"<<endl;
    
    
    vector<vector<vector<double> > > OW_vals_Vec;
    vector<vector<vector<double> > > u_vals_Vec;
    vector<vector<vector<double> > > v_vals_Vec;
    vector<vector<vector<double> > > saltVals_Vec;
    vector<vector<vector<double> > > tempVals_Vec;

    vector<vector<double> > ZetaVals_Vec;
    vector<vector<double> > hVals_Vec;
    vector<vector<double> > lat_rho_Vec;
    vector<vector<double> > lon_rho_Vec;
    vector<vector<double> > pm_Vec;
    vector<vector<double> > pn_Vec;
    vector <double> Cs_r_Vec;
    vector <double> S_rho_Vec;


 
    S_rho_Vec.resize(s_r);
    Cs_r_Vec.resize(s_r);

    lat_rho_Vec.resize(xi_rho);  //  xi_rho = 722
    for (int i = 0; i < xi_rho; ++i)
    {
        lat_rho_Vec[i].resize(eta_rho);// eta_rho = 362   
    }
    
    
    lon_rho_Vec.resize(xi_rho);  //  xi_rho = 722
    for (int i = 0; i < xi_rho; ++i)
    {
        lon_rho_Vec[i].resize(eta_rho);// eta_rho = 362       
    }
    
    
    pm_Vec.resize(xi_rho);  //  xi_rho = 722
    for (int i = 0; i < xi_rho; ++i)
    {
        pm_Vec[i].resize(eta_rho);// eta_rho = 362
    }
      
    
    pn_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        pn_Vec[i].resize(eta_rho);    
    }
    
    
    OW_vals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        OW_vals_Vec[i].resize(eta_rho);
        
        for (int j = 0; j < eta_rho; ++j)
            OW_vals_Vec[i][j].resize(s_r);
    }    
    
    
    u_vals_Vec.resize(xi_rho-1);
    for (int i = 0; i < xi_rho-1; ++i)
    {
        u_vals_Vec[i].resize(eta_rho);
        
        for (int j = 0; j < eta_rho; ++j)
            u_vals_Vec[i][j].resize(s_r);
    }
    
      
    v_vals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        v_vals_Vec[i].resize(eta_rho-1);
        
        for (int j = 0; j < eta_rho-1; ++j)
            v_vals_Vec[i][j].resize(s_r);
    }


    saltVals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        saltVals_Vec[i].resize(eta_rho);
        
        for (int j = 0; j < eta_rho; ++j)
            saltVals_Vec[i][j].resize(s_r);
    }


    tempVals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        tempVals_Vec[i].resize(eta_rho);
        
        for (int j = 0; j < eta_rho; ++j)
            tempVals_Vec[i][j].resize(s_r);
    }
    

    
    hVals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        hVals_Vec[i].resize(eta_rho);
    }
    
    
    ZetaVals_Vec.resize(xi_rho);
    for (int i = 0; i < xi_rho; ++i)
    {
        ZetaVals_Vec[i].resize(eta_rho);
    }
     

    
    long dummyCounter = 0;
    

    for(long y=0; y < s_r ; y++)
    {
        S_rho_Vec[y] = s_rVals[y];
    }
    
    
    for(long y=0; y < s_r ; y++)
    {
        Cs_r_Vec[y] = Cs_r[y];   
    }
    

    dummyCounter = 0;  
    for(long z=0; z < s_r ; z++)
    {
        
        for(long y=0; y < eta_rho ; y++)
        {
            for(long x=0; x < xi_rho ; x++)
            {
                
                tempVals_Vec[x][y][z]= tempVals[dummyCounter];
                dummyCounter++;
            }
        }
    }
       

    dummyCounter = 0;  
    for(long z=0; z < s_r ; z++)
    {
        
        for(long y=0; y < eta_rho ; y++)
        {
            for(long x=0; x < xi_rho ; x++)
            {
                
                saltVals_Vec[x][y][z]= saltVals[dummyCounter];
                dummyCounter++;
            }
        }
    }

    
    dummyCounter = 0;
    for (long z=0; z < s_r ; z++)
    {
        
        for(long y=0; y < eta_rho ; y++)
        {
            for(long x=0; x < xi_rho-1 ; x++)
            {
                if (u_Vals[dummyCounter]< 100000000)
                    u_vals_Vec[x][y][z]= u_Vals[dummyCounter];
                else
                    u_vals_Vec[x][y][z]= 0;
                dummyCounter++;
            }
        }
    }
    

    dummyCounter = 0;
    for(long z=0; z < s_r ; z++)
    {
        
        for(long y=0; y < eta_rho-1 ; y++)
        {
            for(long x=0; x < xi_rho ; x++)
            {
                if (v_Vals[dummyCounter]< 100000000)
                    v_vals_Vec[x][y][z]= v_Vals[dummyCounter];
                else
                    v_vals_Vec[x][y][z]= 0;
                dummyCounter++;
            }
        }
    }  
    
    dummyCounter = 0;
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            
            hVals_Vec[x][y]= h[dummyCounter];
            dummyCounter++;
            //cout<<"hVals_Vec["<<x<<"]["<<y<<"]=["<<hVals_Vec[x][y] <<"]"<<endl;
        }
    }
    
    dummyCounter = 0;
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            if (pm[dummyCounter]< 100000000)
                pm_Vec[x][y]= pm[dummyCounter];
            else
                pm_Vec[x][y]= 0;
            dummyCounter++;
        }
    }
    
    
    dummyCounter = 0;
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            if (pn[dummyCounter]< 100000000)
                pn_Vec[x][y]= pn[dummyCounter];
            else
                pn_Vec[x][y]= 0;
            dummyCounter++;
        }
    }
    


    delete[] u_Vals;
    delete[] v_Vals;
    delete[] h;
    delete[] pm;
    delete[] pn;
    delete[] s_rVals;
    delete[] Cs_r;
    delete[] tempVals;
    delete[] saltVals;
 
    
    
    dummyCounter = 0;
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            if (ZetaVals[dummyCounter] < 100000)
                ZetaVals_Vec[x][y]= ZetaVals[dummyCounter];
            else
                ZetaVals_Vec[x][y]= 0;
            
           // cout<<"ZetaVals_Vec["<<x<<"]["<<y<<"]=["<<ZetaVals_Vec[x][y] <<"]"<<endl;
            
            dummyCounter++;
        }
    }
    
    delete[] ZetaVals;
    
    
    dummyCounter = 0;
    
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            
            lat_rho_Vec[x][y]= lat_rho[dummyCounter];
            dummyCounter++;
            
        }
    }
    
    
    delete[] lat_rho;
    
    
    dummyCounter = 0;
    for(long y=0; y < eta_rho; y++)
    {
        for(long x=0; x < xi_rho ; x++)
        {
            
            lon_rho_Vec[x][y]= lon_rho[dummyCounter];
            dummyCounter++;
        }
    }
  
    
    delete[] lon_rho;
    
  
    
    static int dims[3] = { xi_rho, eta_rho, s_r };
    vtkDataSet * returnVal;
    vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
    sgrid->SetDimensions( dims );
    
    vtkFloatArray *vectors = vtkFloatArray::New();
    int numberofComponents = 1;

    vectors->SetNumberOfComponents(5);
    vectors->SetNumberOfTuples( dims[0] * dims[1] * dims[2] );
    
    vtkPoints *points = vtkPoints::New();
    points->Allocate( dims[0] * dims[1] * dims[2] );
    double dummyMax = 0;
    double dummyMin = 1e+37;
    double currentXVal, currentYVal, currentZVal,z_o;
    double currentS_rho, currentCs_r;
    

    //cout<<"Before the loop!"<<endl;
    long omageIndex = 0;
    double maxZval=-1000;
    double minZval=1000;
    
    double maxYval=-1000;
    double minYval=1000000;
    
    double maxXval=-1000;
    double minXval=100000;

    float OW_total = 0;






    for(long z=0; z < s_r ; z++)
    {
        for(long y=0; y < eta_rho ; y++)
        {
            for(long x=0; x < xi_rho ; x++)
            {
                //    compute z and  values here first.....
                currentS_rho = S_rho_Vec [z];
                currentCs_r = Cs_r_Vec[z];

    
                z_o = ( (hc * currentS_rho ) + ( hVals_Vec[x][y] * currentCs_r ) ) / (hc+hVals_Vec[x][y]);             

                currentZVal = ZetaVals_Vec[x][y] +  (ZetaVals_Vec[x][y] + hVals_Vec[x][y]) * z_o;
                currentZVal = currentZVal / 1000;
               
                if (currentZVal >100000)
                    currentZVal = 10;
                
                if (currentZVal < minZval)
                    minZval = currentZVal;

                
                if ((currentZVal < 10000) && (currentZVal > maxZval))
                    maxZval = currentZVal;

                currentXVal = lon_rho_Vec[x][y];  // lon_rho[dummyCounter];  lon_rho lat_rho
                currentYVal = lat_rho_Vec[x][y];  // lat_rho[dummyCounter];
                
             
                
                if (currentXVal < minXval)
                    minXval = currentXVal;
                
                
                if (currentXVal > maxXval)
                    maxXval = currentXVal;
      
                
                if (currentYVal < minYval)
                    minYval = currentYVal;
                
                
                if (currentYVal > maxYval)
                    maxYval = currentYVal;
                
                
                points->InsertNextPoint(currentXVal, currentYVal, currentZVal);
              
                
                // compute OW and  values here first.....
              
                double S_n, S_s, w_val,  d_x, d_y, d_u_x, d_u_y, d_v_x, d_v_y ;
                
                if(x < xi_rho - 2)
                {
                    d_u_x = u_vals_Vec[ x+1 ][y][z]  -  u_vals_Vec[x][y][z] ;
                    if(y >= eta_rho-1)
                        d_v_x = v_vals_Vec[ x+1 ][y-2][z]  -  v_vals_Vec[x][y-2][z];
                    else
                        d_v_x = v_vals_Vec[ x+1 ][y][z]  -  v_vals_Vec[x][y][z];         
                }
                
                else if(x == xi_rho-1)
                {
                    d_u_x = u_vals_Vec[x-1][y][z] - u_vals_Vec[x-2][y][z];
                    
                    if(y >= eta_rho-1)
                        d_v_x = v_vals_Vec[ x-1 ][y-2][z]  -  v_vals_Vec[x-2][y-2][z];
                    else
                        d_v_x = v_vals_Vec[ x-1 ][y][z]  -  v_vals_Vec[x-2][y][z];                                         
                }
                
                else
                {
                    d_u_x = u_vals_Vec[x-2][y][z] - u_vals_Vec[x-3][y][z];
                    if(y >= eta_rho-1)
                        d_v_x = v_vals_Vec[ x-1 ][y-2][z]  -  v_vals_Vec[x-2][y-2][z];
                    else
                        d_v_x = v_vals_Vec[ x-1 ][y][z]  -  v_vals_Vec[x-2][y][z];
                }
                
                
                
                if(y < eta_rho - 2)
                {
                    d_v_y = v_vals_Vec[ x ][ y+1 ][ z ] - v_vals_Vec[ x ][ y ][ z ];
                    
                    if(x >= xi_rho-1)
                        d_u_y = u_vals_Vec[ x-2 ][y+1][z]  -  u_vals_Vec[x-2][y][z];
                    else
                        d_u_y = u_vals_Vec[ x ][y+1][z]  -  u_vals_Vec[x][y][z];

                    
                }
                else if(y == eta_rho - 1)
                {
                    d_v_y = v_vals_Vec[ x ][ y-1 ][ z ] - v_vals_Vec[ x ][ y-2 ][ z ];                   
                    
                    if(x >= xi_rho-1)
                        d_u_y = u_vals_Vec[ x-2 ][y-1][z]  -  u_vals_Vec[x-2][y-2][z];
                    else
                        d_u_y = u_vals_Vec[ x ][y-1][z]  -  u_vals_Vec[x][y-2][z];
   
                }
                else
                {
                    d_v_y = v_vals_Vec[ x ][ y-2 ][ z ] - v_vals_Vec[ x ][ y-3 ][ z ];
                    if(x >= xi_rho-1)
                        d_u_y = u_vals_Vec[ x-2 ][y-2][z]  -  u_vals_Vec[x-2][y-3][z];
                    else
                        d_u_y = u_vals_Vec[ x ][y-2][z]  -  u_vals_Vec[x][y-3][z];
             
                }
                
                
                
                S_n =  ( d_u_x * pm_Vec[x][y] ) - (d_v_y * pn_Vec[x][y]);
                
                S_s =    (d_v_x * pm_Vec[x][y]) + (d_u_y * pn_Vec[x][y]);
                
                w_val =  (d_v_x * pm_Vec[x][y]) - (d_u_y * pn_Vec[x][y]);
                
                
                OW_vals_Vec[x][y][z] =  ( S_n * S_n ) + ( S_s * S_s ) - ( w_val * w_val );

                OW_total = OW_total + OW_vals_Vec[x][y][z];


                double currentPointVal = OW_vals_Vec[x][y][z];
                
                if (currentPointVal > 1e+36)
                    currentPointVal = 0;
                
                if (currentPointVal > dummyMax )
                    dummyMax = currentPointVal;
                else if (currentPointVal < dummyMin )
                    dummyMin = currentPointVal;
                            

                double currentPointTemp = tempVals_Vec[x][y][z];

                if (currentPointTemp > 100)
                   currentPointTemp = 0;


                double currentPointSalt = saltVals_Vec[x][y][z];

                if (currentPointSalt > 100)
                   currentPointSalt = 0;


                double currentPointU = 0;

                if (x<xi_rho-1)
                   currentPointU = u_vals_Vec[x][y][z]; 
                else
                   currentPointU = u_vals_Vec[x-1][y][z]; 

                if (currentPointU > 200)
                   currentPointU = 0;


                double currentPointV = 0;

                if (y<eta_rho-1)
                   currentPointV = v_vals_Vec[x][y][z]; 
                else
                   currentPointV = v_vals_Vec[x][y-1][z]; 

                if (currentPointV > 200)
                   currentPointV = 0;
                 
                 
                 double currentPointOW = 0;
                 
                if (currentPointVal < -6.31142e-11)
                   currentPointOW = 100;
                   
                else
                	 currentPointOW = 0;


                float * v = new float[5];
                
                v[0] = (float) currentPointOW;
                v[1] = (float) currentPointTemp;
                v[2] = (float) currentPointSalt;
                v[3] = (float) currentPointU;
                v[4] = (float) currentPointV;

                            
                vectors->InsertTuple(omageIndex,v);         
                delete [] v;
                omageIndex++;

 
                if (z==0 && currentXVal>283 && currentYVal>30 && currentYVal<43)
                   fprintf(fpout1,"%f %f %f\n", currentXVal, currentYVal, ZetaVals_Vec[x][y]);             
            }
            
        }       
    }

    fclose(fpout1);



    float OW_mean = OW_total/(s_r*eta_rho*xi_rho);
    float OW_sumstd = 0;


    for(long z=0; z < s_r ; z++)
        for(long y=0; y < eta_rho ; y++)
            for(long x=0; x < xi_rho ; x++)
               OW_sumstd =  OW_sumstd + (OW_vals_Vec[x][y][z] - OW_mean) * (OW_vals_Vec[x][y][z] - OW_mean);
                
    float OW_std = sqrt(OW_sumstd/(s_r*eta_rho*xi_rho));



    *threshold1 = (float) -0.2*OW_std;
    
    cout << (float) -0.2*OW_std << endl;
    
    sgrid->SetPoints(points);
    points->Delete();

    sgrid->GetPointData()->SetScalars(vectors);

    vectors->Delete();
    
    
    // Write file

    string datafile = FileName.substr(0, dotplace) + ".vts";

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =   vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(datafile.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(sgrid);
#else
    writer->SetInputData(sgrid);
#endif
    writer->Write();
   
  
    returnVal = sgrid;
    
    return returnVal; 
}




//-------------------parseConfigFile  FUNCTION----------
//
//------------------------------------------------------

void parseConfigFile(string &base_GeneratedTrackFileName, string &Datapath, string &FileBaseName, string &FileExtention, string Configfilename, int &InitialtimeStep, int &FinaltimeStep, int &SmallestObjVol, int &TimePrecision, int &TimeIncrement, long &x_dim, long &y_dim, long &z_dim, long &x0_dim, long &y0_dim, long &z0_dim,  long &x1_dim, long &y1_dim, long &z1_dim )
{
    // this function reads each individual variable indicated within the Config file and returns them....
    
    string output, datapath, base_filename, initialTimeStep, finalTimeStep,timestepIncrement,fileExtention, timestepPrecision, smallestObjVolToTrack, xdim, ydim, zdim, x0dim, y0dim, z0dim, x1dim, y1dim, z1dim, FileGeneratedTrackFileName;
    cout<<" -----------  Reading the "<<Configfilename<< " Config file ---------"<<endl;
    datapath = string("DATA_FILES_PATH:").c_str();
    base_filename = string("FILE_BASE_NAME:").c_str();
    initialTimeStep = string("INITIAL_TIME_STEP:").c_str();
    finalTimeStep = string("FINAL_TIME_STEP:").c_str();
    timestepIncrement = string("TIME_STEP_INCREMENT:").c_str();
    timestepPrecision = string("TIME_STEP_PRECISION:").c_str();
    smallestObjVolToTrack = string("SMALLEST_OBJECT_VOLUME_TO_TRACK:").c_str();
    fileExtention = string("FILE_EXTENSION:").c_str();
    x1dim =string("Lon_Max:").c_str();
    y1dim =string("Lat_Max:").c_str();
    z1dim =string("Alt_Max:").c_str();
    x0dim =string("Lon_Min:").c_str();
    y0dim =string("Lat_Min:").c_str();
    z0dim =string("Alt_Min:").c_str();
    xdim =string("X_Dim:").c_str();
    ydim =string("Y_Dim:").c_str();
    zdim =string("Z_Dim:").c_str();
    FileGeneratedTrackFileName =string("GENERATED_FILES_PATH:").c_str();
    
    ifstream myreadfile;
    myreadfile.open(Configfilename.c_str(),ios_base::out);
    if (myreadfile.is_open())
    {
        
        while (!myreadfile.eof())
        {
            myreadfile >> output;
            
            if ( output == datapath)
            {
                myreadfile >> output;
                Datapath = output;
            }
            
            if ( output == FileGeneratedTrackFileName)
            {
                myreadfile >> output;
                base_GeneratedTrackFileName = output;
            }
            
            
            if ( output == fileExtention)
            {
                myreadfile >> output;
                FileExtention = output;
            }
            
            if ( output == base_filename)
            {
                myreadfile >> output;
                FileBaseName = output;
            }
            if ( output == initialTimeStep)
            {
                myreadfile >> output;
                InitialtimeStep = atoi(output.c_str());
            }
            if ( output == finalTimeStep)
            {
                myreadfile >> output;
                FinaltimeStep = atoi(output.c_str());
            }
            if ( output == timestepIncrement)
            {
                myreadfile >> output;
                TimeIncrement = atoi(output.c_str());
            }
            
            if ( output == timestepPrecision)
            {
                myreadfile >> output;
                TimePrecision = atoi(output.c_str());
            }
            

            
            if ( output == xdim)
            {
                myreadfile >> output;
                x_dim = atoi(output.c_str());
            }
            
            if ( output == ydim)
            {
                myreadfile >> output;
                y_dim = atoi(output.c_str());
            }
            if ( output == zdim)
            {
                myreadfile >> output;
                z_dim = atoi(output.c_str());
            }
            
            
            if ( output == x1dim)
            {
                myreadfile >> output;
                x1_dim = atoi(output.c_str());
            }
            
            if ( output == y1dim)
            {
                myreadfile >> output;
                y1_dim = atoi(output.c_str());
            }
            if ( output == z1dim)
            {
                myreadfile >> output;
                z1_dim = atoi(output.c_str());
            }
            
            
    
            
            if ( output == x0dim)
            {
                myreadfile >> output;
                x0_dim = atoi(output.c_str());
            }
            
            if ( output == y0dim)
            {
                myreadfile >> output;
                y0_dim = atoi(output.c_str());
            }
            if ( output == z0dim)
            {
                myreadfile >> output;
                z0_dim = atoi(output.c_str());
            }    
            
            if ( output == smallestObjVolToTrack)
            {
                myreadfile >> output;
                SmallestObjVol = atoi(output.c_str());
            }           
        }

        myreadfile.close();
    }

    else
    {
        cout<< "Cannot open the Config.txt File.!!!"<<endl;
    }

    cout<<" Data Path: "<< Datapath <<endl;
    cout<<" Generated File Path: "<< base_GeneratedTrackFileName <<endl;
    cout<<" Base Filename: "<< FileBaseName <<endl;
    cout<<" File Extension: "<< FileExtention<<endl;
    cout<<" Initial TimeStep: "<< InitialtimeStep <<endl;
    cout<<" Final TimeStep: "<< FinaltimeStep <<endl;
    cout<<" Smallest Eddy Volume: "<< SmallestObjVol <<endl;
}





int CreateListFile(string trk_dir, string infile,string &listfile,int InitialtimeStep, int FinaltimeStep,int TimeIncrement,int TimePrecision )
{
    // this function creates a txt file including all the datafilenames that are about to be processed.. (This function is used in the tracking related functions).
    // This function can be removed with the appropriate changes in the future. (A future student can consider this change as a starting project).
    
    int p1,p2;
    p1=infile.rfind('/');
    if (p1==string::npos)
    {
        cout<< "Please use an absolute path when specifying the filename.\n"<< std::endl;
        return 0;
    }
    
    //trk_dir=infile.substr(0,p1+1)+"GENERATED_TRACK_FILES/";
    p2=infile.rfind('.');
    string datafile=infile.substr(p1+1,p2-p1-1); // h5file=bvector
    struct stat st;
    // if trk_dir doesn't exist
    if(stat(trk_dir.c_str(),&st))
    {
        if(mkdir(trk_dir.c_str(),0777)==-1)
        {
            cout << "Unable to create directory: " << trk_dir << std::endl;
            return 0;
        }
        // else if trk_dir exists, but is not a directory
        else if(!(st.st_mode & S_IFDIR))
        {
            cout << "Found a file instead of a directory. Please rename this file, then try again: "<< trk_dir << std::endl;
            return 0;
        }
    }
    // else, use existing directory
    //1.8 create list file and write to list file
    //construct the list file name
    int comp = 0; // ### I have no clue why this comp is used : Rohini
    string listf = trk_dir + datafile + "_" + "comp" + TtoS<int>((int)comp) +
    "_" + TtoS<int>(InitialtimeStep) + "_" + TtoS<int>(FinaltimeStep) +
    "_" + TtoS<int>(TimeIncrement) + ".list";
    //create the list file
    ofstream fp;
    fp.open(listf.c_str(),ios::out|ios::trunc);
    if (!fp.is_open())
    {
        cout << "Unable to open list file: "<< listf << endl;
        return 0;
    }
    //write the directory path
    fp << trk_dir << endl;
    for(register int i=InitialtimeStep;i<FinaltimeStep+1;i+=TimeIncrement)
    {
        string temp;
        temp=datafile+precision_time(i,TimePrecision);
        fp<<temp<<endl;
    }
    fp.close();
    listfile = listf;
    return 1;
}



bool ReadList(char*  ListFileName,vector<string>& time_file, FILETYPE const which, string listfile)
{
    // this function reads the file..
    
    char MainDir[ Consts::MAXLEN];
    char filename[ Consts::MAXLEN];
    FILE *fp;
    if ((fp=fopen(ListFileName, "r"))==NULL)
    {
        fprintf(stderr, "\nCannot open %s file! \n", ListFileName);
        return false;
    }
    
    fscanf(fp, "%s\n", MainDir);
    char buf[ Consts::MAXLEN];
    while (fscanf(fp, "%s\n", filename)!=-1)
    {
        string str(filename);
        
        // if read in a line with only spaces, continue
        std::remove(str.begin(),str.end(),' ');
        if(str.size()==0)
            continue;
        
        // here I need not judge whether the last character of MainDir is '/' or not, because in the spawning version, the listfile is generated by network and I always add '/' behind the string.
        strcpy(filename,str.c_str());
        if(which==_POLY_)
            sprintf(buf, "%s%s.day", MainDir, filename);
        else if(which==_UOCD_)
            sprintf(buf, "%s%s.uocd", MainDir, filename);
        else if(which==_ATTR_)
            sprintf(buf, "%s%s.attr", MainDir, filename);
        else if(which==_BASE_)
            sprintf(buf, "%s%s", MainDir, filename);
        
        time_file.push_back(string(buf));
    }
    fclose(fp);
    return true;
}


bool ReadTrak(string const baseName,Frame& frm)
{
    
    // this file reads the information within the given baseName.
    string const trackFName(baseName+".trak");
    ifstream fp(trackFName.c_str(), ios_base::in);
    if(!fp.is_open())
    {
        cout<<"ReadTrak1:cannot open "<<trackFName<<endl;
        return false;
    }
    string  buf;
    string fName;
    long volume, numObjs(0), numNodes(0);
    float mass, cx, cy, cz,time, l_x, l_y, l_z, u_x,u_y,u_z,minX,minY,minZ,maxX,maxY,maxZ ;
    long vol, Pac_ID;
    
    while(1)
    {
        getline(fp,buf,'\n');
        if(!fp.good())
            break;
        
        ObjectExtends extt;
        istringstream  bufstr(buf,ios_base::out);
        
        bufstr>>mass>>extt.ObjVol>>extt.c_x>>extt.c_y>>extt.c_z>>extt.l_x>>extt.l_y>>extt.l_z>>extt.u_x>>extt.u_y>>extt.u_z>>minX>>minY>>minZ>>maxX>>maxY>>maxZ>>extt.minVal>>extt.maxVal>>Pac_ID;
           
        volume = (long)extt.ObjVol;
        numNodes += volume;
        //numObjs++;
        frm.objVols.push_back(extt);
    }
    fp.close();
    //     frm.objVols.resize(numObjs);
    frm.nodes.resize(numNodes);
    //cout<<"the frm.nodes.size()"<<frm.nodes.size()<<endl;
    return true;
}



bool ReadOct(string baseName,Frame& frm,int step,  vector<vector<TrackObject> > &objs)
{
    // This function reads .uocd files..
    // These files can be used for volume rendering.. (and currently used in the tracking).   
    
    ifstream fp;
    string const octFName(baseName+".uocd");
    fp.open(octFName.c_str(),ios_base::in);
    if(!fp.good())
        cout << "ReadOct:cannot open uocd file to read\n";
    
    float time = 0.0;
    fp>>time;
    int numObjs = 0;
    fp>>numObjs;
    numObjs = frm.objVols.size();
    //cout<<"ReadOct:numObjs="<<numObjs<<endl;
    //cout<<"ReadOct:frm.objVols.size()="<<frm.objVols.size()<<endl;
    assert(numObjs==frm.objVols.size());
    
    //Need not set lengths of objVols and nodes, which has been done in ReadTrak()
    
    float mass = 0.0,cx=0.0, cy=0.0, cz=0.0, Ixx=0.0, Iyy=0.0, Izz=0.0, Ixy=0.0, Iyz=0.0, Izx=0.0;
    unsigned long numNodes = 0, vol = 0, objID = 0;
    register unsigned long i = 0, j = 0;
    unsigned long x = 0, y = 0 , z = 0,  xsize = 0, ysize = 0, zsize = 0, vertID = 0;
    float xPos = 0, yPos = 0, zPos = 0, val = 0, x1 = 0, y1 = 0, z1 = 0;
    unsigned long k(0);
    
    
    
    for (i=0; i<numObjs; i++)
    {
        fp>>objID>>vol>>mass>>cx>>cy>>cz>>Ixx>>Iyy>>Izz>>Ixy>>Iyz>>Izx \
        >>vertID>>x1>>y1>>z1>>val;
        

        //NOTE: x1,y1,z1 are not local maximum actually! need rectifying
        
        objs[step][objID].existent =  Consts::YES;
        objs[step][objID].time = time;
        objs[step][objID].mass = mass;
        objs[step][objID].volume = vol;
        objs[step][objID].xbar = cx;
        objs[step][objID].ybar = cy;
        objs[step][objID].zbar = cz;
        objs[step][objID].Ixx=Ixx;
        objs[step][objID].Iyy=Iyy;
        objs[step][objID].Izz=Izz;
        objs[step][objID].Ixy=Ixy;
        objs[step][objID].Iyz=Iyz;
        objs[step][objID].Izx=Izx;
        objs[step][objID].maximum=val;
        objs[step][objID].max_x=(int)x1;
        objs[step][objID].max_y=(int)y1;
        objs[step][objID].max_z=(int)z1;
        
        frm.nodes[k].NodeID = vertID;
        frm.nodes[k].ObjID=i;
        k++;
        
        for (j=0; j<vol-1; j++)
        {
            fp>>vertID>>x1>>y1>>z1>>val;
            frm.nodes[k].NodeID = vertID;
            frm.nodes[k].ObjID=i;
            k++;
        }
    }
    fp.close();
    return true;
}




vector<string>
Tokenize(const string &sentence)
{
    // this function reads string and tokenize the items in the string seperated by space..
    
    stringstream ss(sentence);
	istream_iterator<string> it(ss);
	istream_iterator<string> end;
	vector<string> tokens(it,end);
	return tokens;
}


void readFile(string fullFilename, vector<string> &fullFrameInfo)
{
    // this function reads the entire file content into a vector where each line of the file is a string variable
    fullFrameInfo.clear();
    ifstream input1;
    input1.open(fullFilename.c_str());
    if(!input1)
    {
        cout<<"\n readFile: Cannot open "<<fullFilename<<" file "<<endl;
        return;
    }
    string str1;
    getline(input1, str1);
    while(!input1.eof())
    {
        fullFrameInfo.push_back(str1);
        getline(input1, str1);
    }
    input1.close();
    
}



void findObjIDinPreviousFrame(vector<long> & previousObjects,vector<string> &trakTable, long currentObj )
{
    // takes the currentObj number (starts from 1) and checks it in the tracktable file to find it in the second half
    // then returns the objectIDs in the first half in the previousObjects vector..
    // currentObj starts from 1 since it is in traktable...
    // it returns -1 (if the object does not exist in the previous frame) or returns all the possible objectIDs in the previousObjects
    //
    long linenumber(0),k(0), previousObjnum(0),secondframe;
    vector<string> tokensforTrakTable;
    bool found = 0;
    previousObjects.clear();
    
    while((linenumber < trakTable.size()) && !found ) //find the currentObj within the trakTable
    {
        tokensforTrakTable.clear();
        tokensforTrakTable = Tokenize(trakTable[linenumber].c_str()); //get the current line of the trakTable
        k = 0;
        previousObjnum = -1;
        secondframe = 0;
        while ( (!found)  && ( k < tokensforTrakTable.size() ) ) // check along the each line within th traktable for the current object
        {
            // k is the tokensforTrakTable index
            // remember objects start from 0 within the Packets, however they start from 1 within the .trakTable file, consider that!!!...
            if (atoi(tokensforTrakTable[k].c_str()) == currentObj) // we found the currentObj
            {
                for (long secondframe1 = 0; secondframe1 < tokensforTrakTable.size(); secondframe1++)
                {
                    if ((tokensforTrakTable[secondframe1] == "-1") && (secondframe1 < k) ) // check if the found object within the first frame or within the second frame
                    {
                        found = 1;
                        if (secondframe1 != 0 )
                        {
                            int dum1 = secondframe1-1;
                            while ( 0 <= dum1)
                            {
                                previousObjnum = atoi(tokensforTrakTable[dum1].c_str());
                                previousObjects.push_back(previousObjnum);
                                dum1--;
                            }
                        }
                        else
                        {
                            previousObjnum = -1;
                            previousObjects.push_back(previousObjnum);
                        }
                    }
                }
            }
            k++;
        }
        linenumber++;
    }
}





void
ReadCurrentFrameTrackInfo(string trakTableFile, int cycle, vector<string> &trakTable)
{
    
    //cout<<" INSIDE READCURRENTFRAMETRACKINFO !!! the cyle value is:[" <<  cycle <<"]"<<endl;
    trakTable.clear(); // make sure the trakTable vector is clear ......!!
	ifstream input1;
	input1.open(trakTableFile.c_str());
    
	if(!input1)
	{
		cout<<"\n Cannot open "<<trakTableFile<<" file "<<endl;
		return;
	}
    
    string find_frame = "Frame #" +TtoS<int>(cycle);
	string str1;
    
    //cout<<"find frame = "<<find_frame<<endl;
    
	bool frame_found = false;
	
	while(!input1.eof())
	{
		getline(input1, str1);
		
		if(str1[0] == 'F')
		{
			if (str1 == find_frame)
			{
				frame_found = true;
			}
			else
			{
				frame_found = false;
			}
		}
        
		if(frame_found && str1[0] != 'F')
		{
			trakTable.push_back(str1);
		}
	}
}








int BeginFeatureTrack(string base_GeneratedTrackFileName, string currenttimevalue, char *mypolyfile,int &finished, string listfile, int InitialtimeStep, int timestepincrement, int currentTime, int TimePrecision)
{
    vector<string>   time_polyfile; // seq_wr+no_init
    vector<string>   time_filename; // seq_wr+no_init

    char *tempstr;
    
  
    if(!ReadList((char*)listfile.c_str(), time_filename, _BASE_, listfile)) // <<<<<<
    {
        //OutPolyFilename.set_str_val(PolyFilename_.c_str());

        cout << "fail to read list file!" << endl;
        
        return 1;
    }

    if(!ReadList((char*)listfile.c_str(), time_polyfile, _POLY_, listfile)) // <<<<<<
    {

        cout << "fail to read *.day file!" << endl;
        
        return 1;
    }
   


    int i(0);
    int TrackedFr(-1);
    for(vector<string>::const_iterator it=time_polyfile.begin(); it!=time_polyfile.end(); ++it,++i)
    {
        if(mypolyfile== *it)
        {
            TrackedFr=i;
            break;
        }
    }

    string label(time_filename.at(0)); //this is the (full) path to GENERATED_TRACK_FILES folder..
    string label1(time_filename.at(TrackedFr)); //this is the current file's number ...
    string currentfile(mypolyfile);

    currentfile = currentfile.substr(0, currentfile.rfind('/')+1);

    
    int prev_temp_index = ((int)InitialtimeStep+(currentTime-1)*timestepincrement);
    string prevtimevalue = precision_time(prev_temp_index,TimePrecision);
    
    string colormapfile = currentfile + "ColorMap" + currenttimevalue + ".txt";
    string previouscolormapfile = currentfile + "ColorMap" + prevtimevalue + ".txt";

    if(TrackedFr==-1)
        cout  <<"trk:bad list file!\n";

    if(TrackedFr==0 && time_polyfile.size()>1)
    {
        cout<<"this is the first frame, return without tracking\n";
        finished=0;

        
        return 1;
    }

    if(TrackedFr==0 && time_polyfile.size()==1)
    {
        cout<<"only one frame,tracking finished\n";
        finished=1;	 
        return 1;
    }

    if(TrackedFr>1)
    {
        fstream outfile;
        outfile.open(string(label+".trakTable").c_str(),ios_base::in);
        if(!outfile.is_open())
            cout << "please try from the first frame, return...\n";
    }

    Frame t1, t2;

    
    if (!ReadTrak(time_filename.at(TrackedFr-1), t1))      // because the minimum of TRackedFr is 1(if 0, then program should have returned earlier)
        return 0;
    vector<vector<TrackObject> >
    objs( Consts::DEFAULT_TIMESTEP_NUM, vector<TrackObject>( Consts::MAXOBJS));
    
    if (!ReadOct(time_filename.at(TrackedFr-1), t1, TrackedFr-1,objs)) 
        return 0;
    
    sort(t1.nodes.begin(),t1.nodes.end());   
    bound_check<vector<TrackObject> >(objs,TrackedFr-1, Consts::DEFAULT_TIMESTEP_NUM);
    objs.at(TrackedFr-1).at(0).time=(float)TrackedFr-1;
    objs.at(TrackedFr).at(0).time=(float)TrackedFr;
    
    if(!ReadTrak(time_filename.at(TrackedFr), t2))
        return 0;
    if(!ReadOct(time_filename.at(TrackedFr), t2, TrackedFr,objs))
        return 0;
    sort(t2.nodes.begin(),t2.nodes.end());
    

    TrackObjects(label,TrackedFr, (InitialtimeStep + timestepincrement*currentTime), t1, t2, time_polyfile, objs); //<<<<<<<<<
    string traktablefile = label + ".trakTable";
    vector<string> trakTable;    
    ReadCurrentFrameTrackInfo(traktablefile, InitialtimeStep + timestepincrement*currentTime, trakTable);       
    
    //cout<<"=====current TrackedFr & t2.objVols.size()===== "<<TrackedFr<<" & "<<t2.objVols.size()<<endl;
    
    
    if (TrackedFr==time_polyfile.size()-1)
       finished=1;
   
    else
       finished=0;
       
    objs.clear();
    trakTable.clear();
    
    return 1;
}




//-------------------MAIN FUNCTION----------------------
//
//------------------------------------------------------

int main(void)
{
    unsigned long ncells, nnodes;
    int cellpoints, nncomp, count, count1(0), count2(0), found, InitialtimeStep, FinaltimeStep, SmallestObjVol, TimePrecision, TimeIncrement, currentTime;
    long x_dim(0), y_dim(0), z_dim(0), lon_min(0), lat_min(0), alt_min(0), lon_max(0), lat_max(0), alt_max(0);
    int nspace = 3;    //data space dimension
    int numberOfPointsInCell = 8;    // assuming that the cell is cube.. same variable as cellpoints
    string file_name, datapath, FileBaseName, fileextension, listfile, base_output, zetafile;
    unsigned long *node_conn_array = NULL;
    float *coord_array = NULL;
    float *node_data = NULL;
    file_name = string("Config.txt").c_str();
    cellpoints =8;
    
 
    parseConfigFile(base_output, datapath, FileBaseName, fileextension, file_name, InitialtimeStep, FinaltimeStep, SmallestObjVol, TimePrecision, TimeIncrement, x_dim, y_dim, z_dim, lon_min, lat_min, alt_min, lon_max, lat_max, alt_max);
    
    time_t begin1, end1;
    string polyext = ".day";
    string nonamecharacter ="+";
    if (strcmp (FileBaseName.c_str(),nonamecharacter.c_str()) == 0 )
        int a = CreateListFile(base_output,datapath,listfile,InitialtimeStep,FinaltimeStep,TimeIncrement,TimePrecision );
    else
        int a = CreateListFile(base_output, datapath+FileBaseName, listfile, InitialtimeStep, FinaltimeStep, TimeIncrement, TimePrecision );
    
    
    int maxCycleNumber = ((FinaltimeStep-InitialtimeStep)/TimeIncrement) + 1;
    
    for (int currentTime = 0; currentTime <maxCycleNumber; currentTime++ )
    {
        clock_t begin = clock();
        
        string base_filename;
        string trakfile1;
        string attributeFile;
    
        string OutputOcdfile;
        int temp_index = ((int)InitialtimeStep+currentTime*TimeIncrement);
        string currenttimevalue = precision_time(temp_index,TimePrecision);
        if (strcmp (FileBaseName.c_str(),nonamecharacter.c_str()) == 0 )
            base_filename = datapath+currenttimevalue;
        else
            base_filename = datapath+FileBaseName+currenttimevalue;
        
        file_name = base_filename + fileextension;
        
        string base_GeneratedTrackFileName = base_output; //datapath + "GENERATED_TRACK_FILES/";
        
        string mypolyfile;
        if (strcmp (FileBaseName.c_str(),nonamecharacter.c_str()) == 0 )
        {
            mypolyfile = base_GeneratedTrackFileName + currenttimevalue + polyext;
            trakfile1 = base_GeneratedTrackFileName  + currenttimevalue + ".trak";
            attributeFile = base_GeneratedTrackFileName+currenttimevalue + ".attr";     
            OutputOcdfile = base_GeneratedTrackFileName + currenttimevalue + ".uocd";
        }

        else
        {
            trakfile1 = base_GeneratedTrackFileName  + FileBaseName +  currenttimevalue + ".trak";
            mypolyfile = base_GeneratedTrackFileName + FileBaseName +  currenttimevalue + polyext;
            attributeFile = base_GeneratedTrackFileName + FileBaseName + currenttimevalue + ".attr";       
            OutputOcdfile = base_GeneratedTrackFileName + FileBaseName +currenttimevalue + ".uocd";
            base_GeneratedTrackFileName = base_GeneratedTrackFileName + FileBaseName;
        }
         
        float thresh1;
        vtkDataSet * in_ds = ReadNcDataFile(base_output, file_name, x_dim, y_dim, z_dim, &thresh1);       
        nnodes = x_dim * y_dim * z_dim;
        
        int doVolRender = 1;
        
        cout<< endl;
        cout<< endl;
        cout<< endl;
        cout<<"-------- Computation of the "<< currentTime+1 << "th file begins now ! ---------------"<<endl;
        
        vtkCellTypes *types = vtkCellTypes::New();
        in_ds->GetCellTypes(types);
        int cell_types = types->GetNumberOfTypes();
        if(cell_types > 1)
        {
            cout << " warning: There are multiple cell types! Total number of Cell types is: " << cell_types << "\n only 1st type will be processed \n";
        }
        types->Delete();
        int cell_type = in_ds->GetCellType(0);
        
        
        if (cell_type == VTK_VOXEL)
           cout << "Cell Type: voxel" << endl;

        else if (cell_type == VTK_QUAD)
           cout << "Cell Type: quad" << endl;

        else if (cell_type == VTK_HEXAHEDRON)
           cout << "Cell Type: hexahedron" << endl;

        else if (cell_type == VTK_PIXEL)
           cout << "Cell Type: pixel" << endl;

        else
        {
            cout << "Only voxel, quad or hexahedron type cell can be processed, cell_type is:["<< cell_type<< "]" <<endl;
            return 0;
        }

        vtkDataSet *outDS = NULL;

        vtkIdType ncells = in_ds->GetNumberOfCells();
       
        
        bool *ProcessedCell=0;
        ProcessedCell = new bool[ncells];
        for (long i=0; i<ncells ; i++ )
            ProcessedCell[i] =false;

        int ncomp = in_ds->GetPointData()->GetNumberOfComponents(); //this is the number vector /scalars etc
        vtkDataArray *vtk_node_data = NULL;

        //cout << "Number of Components:" << ncomp << endl;
       
        vtk_node_data = in_ds->GetPointData()->GetScalars();
       
  
        std::set<vtkIdType>::iterator it1;
        std::set<vtkIdType>::iterator it2;
        int aa = 8;              // (int) cellPointIds->GetNumberOfIds(); // this is 8 for a cube /rectilinear data -- (=numberofcellpoints).
        long firstobj = -1;

        
        for (unsigned long i  = 0; i < ncells; i++) // go for each cell
        {
            if (ProcessedCell[i] == false)
            {
                //ProcessedCell[i] = true;
                
                vtkIdType cellId = (vtkIdType) i;  //<-- this is a cell's id. not a point's Id.
                
                
                std::set<vtkIdType> ThisObjectsMemberPointIds;
                std::set<vtkIdType> CandidateNeigborCellIdsToProcess;
                std::set<vtkIdType> ThisObjectsCellIds;
     
                CandidateNeigborCellIdsToProcess.clear();
                ThisObjectsMemberPointIds.clear();
                ThisObjectsCellIds.clear();
         
                CandidateNeigborCellIdsToProcess.insert(cellId);
                         
                
                double MaxXDimValue = -999999999999999;
                double MaxYDimValue = -999999999999999;
                double MaxZDimValue = -999999999999999;
                
                double MinXDimValue = 999999999999999;
                double MinYDimValue = 999999999999999;
                double MinZDimValue = 999999999999999;
                
                
                while(CandidateNeigborCellIdsToProcess.size()>0)
                {
                    it1 = CandidateNeigborCellIdsToProcess.begin();
                    cellId = *it1; 
                    CandidateNeigborCellIdsToProcess.erase (cellId);
                    if(ProcessedCell[(unsigned long)cellId] == false)
                    {
                        ProcessedCell[(unsigned long)cellId] = true;
                        vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
                        in_ds->GetCellPoints(cellId, cellPointIds); //<--- this one gets the points of the current cell (cellId)
                        for(vtkIdType j = 0; j < aa ; j++) // for each cell point in this cell
                        {
                            vtkIdType pointindex = cellPointIds->GetId(j);
                            double* temp_val = vtk_node_data->GetTuple(pointindex);
                              
                            if (temp_val[0]<=thresh1) 
                            {        
                                double *temp_coord3D = new double[3];
                                in_ds->GetPoint(pointindex,temp_coord3D);
                                
                                if (temp_coord3D[0] > MaxXDimValue)
                                    MaxXDimValue = temp_coord3D[0];
                                
                                if (temp_coord3D[1] > MaxYDimValue)
                                    MaxYDimValue = temp_coord3D[1];
                                
                                if (temp_coord3D[2] > MaxZDimValue)
                                    MaxZDimValue = temp_coord3D[2];
                                                               
                                if (temp_coord3D[0] < MinXDimValue)
                                    MinXDimValue = temp_coord3D[0];                                                        
                                
                                if (temp_coord3D[1] < MinYDimValue)
                                    MinYDimValue = temp_coord3D[1];                                
                                
                                if (temp_coord3D[2] < MinZDimValue)
                                    MinZDimValue = temp_coord3D[2];

                                delete[] temp_coord3D;
                     

                                ThisObjectsCellIds.insert(cellId);
                                ThisObjectsMemberPointIds.insert(pointindex);
                                vtkSmartPointer<vtkIdList> thispointsCellIds = vtkSmartPointer<vtkIdList>::New();
                                in_ds-> GetPointCells(pointindex, thispointsCellIds);

                                for (vtkIdType k = 0; k < thispointsCellIds->GetNumberOfIds(); k++)                              
                                    CandidateNeigborCellIdsToProcess.insert(thispointsCellIds->GetId(k));                             
                            }
                            
                        }
                        
                    }
                    
                }

          
                
                bool conditionforThisDataType = ( (ThisObjectsMemberPointIds.size()>=SmallestObjVol) && (MinXDimValue > lon_min) && (MaxXDimValue < lon_max) && (MinYDimValue > lat_min) && (MaxYDimValue < lat_max) && (MinZDimValue > alt_min) && (MinZDimValue < -0.05) && (MaxZDimValue > -0.02) );

                                            
                if (conditionforThisDataType)  //this is the extracted object!!!...
                {
                    firstobj++;
                    unsigned long thisObjectsVolume = ThisObjectsMemberPointIds.size();                   
                           
                    vtkIdType *ObjectCellsArray = new vtkIdType[ThisObjectsCellIds.size()];
                    //vtkIdType *CellPointIdsArray = new vtkIdType[thisObjectsVolume];
                    long tt=0;
                    
                                      
                    // convert from set to an array..
                    tt=0;
                    for (std::set<vtkIdType>::iterator it=ThisObjectsCellIds.begin(); it!=ThisObjectsCellIds.end(); ++it)
                    {
                        ObjectCellsArray[tt] = *it;
                        tt++;
                    }
                    ThisObjectsCellIds.clear();

                    int dummycounter = 0;
                    int founditem = 0;
                    unsigned long totalnumberofcellsinthisobj = tt;//ThisObjectsCellIds.size();

         
                    //create a new unstructured data set for this object only..
                    
                    std::set<vtkIdType> ThisObjectsCompleteCellPointIds;
                    ThisObjectsCompleteCellPointIds.clear();
                    
                    
                    vtkSmartPointer<vtkIdList> cellPointIds1 = vtkSmartPointer<vtkIdList>::New();
                    for(unsigned long hm = 0; hm < totalnumberofcellsinthisobj; hm++)
                    {
                        in_ds->GetCellPoints(ObjectCellsArray[hm], cellPointIds1);
                        for(vtkIdType j = 0; j < aa ; j++)
                        {
                            vtkIdType pointindex1 = cellPointIds1->GetId(j);
                            ThisObjectsCompleteCellPointIds.insert(pointindex1);
                        }
                    }
                    
                    unsigned long cursor=0;
                    cursor = ThisObjectsCompleteCellPointIds.size();
                    vector<vtkIdType> PointIdsArray1;
                    vector<vtkIdType> PointIdsArray;
                    PointIdsArray1.clear();
                    PointIdsArray.clear();
                    PointIdsArray1.resize(cursor);
                    PointIdsArray.resize(thisObjectsVolume);
                
                    long tt1=0;
                    for (std::set<vtkIdType>::iterator it=ThisObjectsCompleteCellPointIds.begin(); it!=ThisObjectsCompleteCellPointIds.end(); ++it)
                    {
                        PointIdsArray1[tt1] = *it;
                        tt1++;
                    }
                    
                    
                    long tt2=0;
                    for (std::set<vtkIdType>::iterator it=ThisObjectsMemberPointIds.begin(); it!=ThisObjectsMemberPointIds.end(); ++it)
                    {
                        PointIdsArray[tt2] = *it;
                        tt2++;
                    }
                    
      
                    float *coords=0; //convert this back to double soon. :)
                    double *nodevals = new double[cursor];
                    double *nodevals1 = new double[cursor];
                    double *nodevals2 = new double[cursor];
                    double *nodevals3 = new double[cursor];
                    double *nodevals4 = new double[cursor];
                    unsigned long *connects=0;
                    connects=new unsigned long[totalnumberofcellsinthisobj*cellpoints];
                    double *temp_coord = new double[3];
                    coords=new float[cursor*nspace];

                    unsigned long j1=0;
                    unsigned long j2=0;
                    unsigned long j3=0;
                    double* temp_val;
                    double mass = 0;
                    double Cx=0, Cy=0, Cz=0;
                    double LowerLeft_Corner_x, LowerLeft_Corner_y, LowerLeft_Corner_z;
                    double UpperRight_Corner_x, UpperRight_Corner_y, UpperRight_Corner_z;
                    LowerLeft_Corner_x = 99999999999.99; // in the future assign these values to some constant variable and use that varible here ..
                    LowerLeft_Corner_y = 99999999999.99;
                    LowerLeft_Corner_z = 99999999999.99;
                    UpperRight_Corner_x = -999999999.99;
                    UpperRight_Corner_y = -999999999.99;
                    UpperRight_Corner_z = -999999999.99;
                    double MinVal = 999999999, MaxVal = -999999999;
                    double MinValX = 0, MinValY = 0, MinValZ = 0, MaxValX = 0, MaxValY = 0, MaxValZ = 0;   
                    double MinTemp = 999, MaxTemp = -999, MinSalt = 999, MaxSalt = -999, MinSpeed = 999, MaxSpeed = -999;  
                    
                    vtkIdType maxValNodeID=0;
                    vtkIdType minValNodeID=0;                                        

                    
                    for (vtkIdType kt1=0; kt1< cursor ; kt1++)
                    {
                        vtkIdType thisCurrentPoint = PointIdsArray1[kt1];
                        temp_val = vtk_node_data->GetTuple(PointIdsArray1[kt1]);
                        nodevals[j2] = temp_val[0];
                        nodevals1[j2] = temp_val[1];
                        nodevals2[j2] = temp_val[2];
                        nodevals3[j2] = temp_val[3];
                        nodevals4[j2++] = temp_val[4];                                          
                        
                        in_ds->GetPoint(thisCurrentPoint,temp_coord);
                        
                        
                        double speedsq = temp_val[3]*temp_val[3] + temp_val[4]*temp_val[4];
                        
                        if (temp_val[1] < MinTemp)
                           MinTemp = temp_val[1];
                           
                        if (temp_val[1] > MaxTemp)
                           MaxTemp = temp_val[1];
                           
                        if (temp_val[2] < MinSalt)
                           MinSalt = temp_val[2];
                           
                        if (temp_val[2] > MaxSalt)
                           MaxSalt = temp_val[2];
                           
                        if (speedsq < MinSpeed)
                           MinSpeed = speedsq;
                           
                        if (speedsq > MaxSpeed)
                           MaxSpeed = speedsq;
                        
                        
                        for (int dim1=0; dim1<3; dim1++)
                            coords[j1++] = temp_coord[dim1];
                       
                        int ismemberofobject = ThisObjectsMemberPointIds.count(thisCurrentPoint);
                        
                        if (ismemberofobject!=0)
                        {
                            Cx += (temp_val[0]*temp_coord[0]);
                            Cy += (temp_val[0]*temp_coord[1]);
                            Cz += (temp_val[0]*temp_coord[2]);
                            mass = mass + temp_val[0];
                            
                            if ( temp_coord[0] < LowerLeft_Corner_x )
                                LowerLeft_Corner_x = temp_coord[0];
                            else if ( temp_coord[0]  > UpperRight_Corner_x )
                                UpperRight_Corner_x = temp_coord[0];
                            
                            if ( temp_coord[1] < LowerLeft_Corner_y )
                                LowerLeft_Corner_y = temp_coord[1];
                            else if ( temp_coord[1]  > UpperRight_Corner_y )
                                UpperRight_Corner_y = temp_coord[1];
                            
                            
                            if ( temp_coord[2] < LowerLeft_Corner_z )
                                LowerLeft_Corner_z = temp_coord[2];
                            else if ( temp_coord[2]  > UpperRight_Corner_z )
                                UpperRight_Corner_z = temp_coord[2];
                            
                            if (MinVal > temp_val[0])
                            {
                                minValNodeID = thisCurrentPoint;
                                MinVal = temp_val[0];
                                MinValX = temp_coord[0];
                                MinValY = temp_coord[1];
                                MinValZ = temp_coord[2];
                            }
                            else if (MaxVal < temp_val[0])
                            {
                                maxValNodeID = thisCurrentPoint;
                                MaxVal = temp_val[0];
                                MaxValX = temp_coord[0];
                                MaxValY = temp_coord[1];
                                MaxValZ = temp_coord[2];
                                
                            }
                        }
                                               
                    }
                    
                    Cx = Cx /mass;
                    Cy = Cy /mass;
                    Cz = Cz /mass;
                   
                    
                    //cout << "Liu: " << MinTemp << ", " << MaxTemp << endl;
            
                    
                    //compute the Second Moments.....        
                    double Ixx=0.0, Iyy=0.0, Izz=0.0, Iyz=0.0, Ixy=0.0, Izx=0.0;
                    double mx = 0.0, my = 0.0, mz = 0.0;
                    
                    for (long k5 = 0; k5< ThisObjectsMemberPointIds.size(); k5++)
                    {
                        vtkIdType thisCurrentPoint = PointIdsArray[k5];
                        temp_val = vtk_node_data->GetTuple(thisCurrentPoint);
                        in_ds->GetPoint(thisCurrentPoint,temp_coord); //temp_coord[0]=x, [1]=y, [2]=z,
                        mx = temp_coord[0]-Cx;
                        my = temp_coord[1]-Cy;
                        mz = temp_coord[2]-Cz;
                        
                        Ixx += temp_val[0]*mx*mx;
                        Iyy += temp_val[0]*my*my; //(val*pow(my,2)/((float)(objPtr->mass)));
                        Izz += temp_val[0]*mz*mz; //(val*pow(mz,2)/((float)(objPtr->mass)));
                        Iyz += temp_val[0]*my*mz; //(val*my*mz/((float)(objPtr->mass)));
                        Ixy += temp_val[0]*mx*my; //(val*mx*my/((float)(objPtr->mass)));
                        Izx += temp_val[0]*mz*mx; //(val*mz*mx/((float)(objPtr->mass)));                       
                    }
                    
                    Ixx = Ixx/(mass);
                    Iyy = Iyy/(mass);
                    Izz = Izz/(mass);
                    Iyz = Iyz/(mass);
                    Ixy = Ixy/(mass);
                    Izx = Izx/(mass);
                    
                    
                    // set the cellpoint ids to a new number starting from 0. :)
                    
                    vtkSmartPointer<vtkIdList> cellPointIds2 = vtkSmartPointer<vtkIdList>::New();
 
                    for (unsigned long i = 0; i < totalnumberofcellsinthisobj;i++)
                    {
                        in_ds->GetCellPoints(ObjectCellsArray[i], cellPointIds2);
                        for(int j = 0; j < cellpoints;j++)
                        {
                            vtkIdType j4 = 0;
                            vtkIdType realpointindex1 = cellPointIds2->GetId(j);
                            int notfound =1;
                            while (notfound)
                            {
                                 if (PointIdsArray1[j4] == realpointindex1)
                                    notfound  = 0;
       
                                 j4++;                           
                            }
                            
                            connects[j3++] = j4-1;
                        }
                    }
                                        
                    delete[] ObjectCellsArray;
                    
                    
                    /********************** Output track file ***********************/
                    /*
                    ofstream Trakfile;
                    
                    long volume = thisObjectsVolume;
                    
                    if (firstobj ==0)
                       Trakfile.open(trakfile1.c_str(),ofstream::out | ofstream::trunc);
                 
                    else           
                       Trakfile.open(trakfile1.c_str(),ofstream::out | ofstream::app);
                        
                  
                    Trakfile << mass << "   "<<volume<< "   " << Cx <<  "   " << Cy << "   " << Cz << "   "<< LowerLeft_Corner_x <<"   ";
                    Trakfile << LowerLeft_Corner_y << "   "   << LowerLeft_Corner_z << "   " << UpperRight_Corner_x << "   " ;
                    Trakfile << UpperRight_Corner_y << "   " <<UpperRight_Corner_z << "   "  << MinValX    <<"   "<< MinValY << "   " << MinValZ << "   " << MaxValX << "   " << MaxValY << "   " << MaxValZ << "   " << MinVal << "   " << MaxVal << endl;
                    Trakfile.close();
                    */


                    long volume = thisObjectsVolume;

                    zetafile = base_GeneratedTrackFileName + currenttimevalue + ".zeta";
                    float spin_orientation = search_zeta(zetafile, Cx, Cy, LowerLeft_Corner_x, LowerLeft_Corner_y, UpperRight_Corner_x, UpperRight_Corner_y);
                    
                    //remove(zetafile.c_str());

                    /********************** OutputAttribute file ***********************/
                    
                    FILE *fpout;
                    //int rc;
                    //char buffer[256];
                    char outAttr[256];
                    strcpy(outAttr,attributeFile.c_str());                 
                    
                    if (firstobj ==0)
                    {
                       if ((fpout = fopen(outAttr, "w"))==NULL)
                          cout << "cannot open outAttr file to write\n";
                    }
                    
                    else
                    {
                       if ((fpout = fopen(outAttr, "a"))==NULL)
                          cout << "cannot open outAttr file to write\n";
                        
                    }
                    
    
                    fprintf(fpout,"------------------------------------------------\n");
                    fprintf(fpout,"object %ld attributes:\n", firstobj);
                    fprintf(fpout,"Max position: (%f, %f, %f) with value: %f\n", MaxValX, MaxValY, MaxValZ, MaxVal);

                    fprintf(fpout,"Node #: %lld\n",maxValNodeID);
                    fprintf(fpout,"Min position: (%f, %f, %f) with value: %f\n", MinValX, MinValY, MinValZ, MinVal);
                    fprintf(fpout,"Node #: %lld\n",minValNodeID);
                    fprintf(fpout,"Integrated content: %f\n",mass);
                    fprintf(fpout,"Sum of squared content values: %f\n",mass*mass); // this is not the real mass sequare definition. Fix this later.
                    fprintf(fpout,"Volume: %ld\n",volume);
                    fprintf(fpout,"Centroid: (%f, %f, %f)\n",Cx, Cy, Cz);
                    fprintf(fpout,"Moment: Ixx = %f\nIyy = %f\nIzz = %f\n", Ixx, Iyy, Izz);
                    fprintf(fpout,"Ixy = %f\nIyz = %f\nIzx = %f\n", Ixy, Iyz, Izx);
                    
                    fclose(fpout);
                    
                    
                    
                    /********************** Output OutputOcd file ***********************/
               
                    FILE *fpout1;
                    char OutOcd[256];
                    strcpy(OutOcd,OutputOcdfile.c_str());
        
                    if(firstobj ==0)
                    {
                        //cout<<"inside if. firstobj= " <<firstobj <<endl;
                        if ((fpout1 = fopen(OutOcd, "w"))==NULL)
                            cout << "cannot open outAttr file to write\n";
                        fprintf(fpout1,"%f\n",(float) currentTime);
                        fprintf(fpout1,"%ld\n", firstobj);//this should be total number of object... we do not really need this. And adding this to the file requires a full obj loop. Not necessary... look into this.
                    }
                    else
                    {
                        if((fpout1 = fopen(OutOcd, "a"))==NULL)
                            cout << "cannot open outAttr file to write\n";
                        
                    }
                    
                    fprintf(fpout1,"%ld\n", firstobj);
                    fprintf(fpout1,"%ld %f %f %f %f\n", thisObjectsVolume, (float) mass, (float) Cx, (float) Cy, (float) Cz );
                    fprintf(fpout1,"%f %f %f %f %f %f\n", (float) Ixx, (float) Iyy, (float) Izz, (float) Ixy, (float) Iyz, (float) Izx);
                    
                    for(long k5 = 0; k5< ThisObjectsMemberPointIds.size(); k5++)
                    {
                        vtkIdType thisCurrentPoint = PointIdsArray[k5];
                        temp_val = vtk_node_data->GetTuple(thisCurrentPoint);
                        in_ds->GetPoint(thisCurrentPoint,temp_coord); //temp_coord[0]=x, [1]=y, [2]=z,
                        

                        fprintf(fpout1,"%6lld %9.6f %9.6f %9.6f  %f\n", thisCurrentPoint, (float) temp_coord[0], (float) temp_coord[1], (float) temp_coord[2], (float) temp_val[0]);
                                               
                    }
                    
                    fclose(fpout1);
                         
                    
                    
                    //------isosurface related ----------------------------------
                    
                    //vtkUnstructuredGrid *out_ds1 = vtkUnstructuredGrid::New();
                    
                    vtkSmartPointer<vtkUnstructuredGrid> out_ds1 = vtkSmartPointer<vtkUnstructuredGrid>::New();
                    //--------Below is Setoutfield -----
                    vtkSmartPointer<vtkFloatArray> pcoords = vtkSmartPointer<vtkFloatArray>::New();
                    pcoords->SetNumberOfComponents(nspace);
                    pcoords->SetNumberOfTuples(cursor);
                    float *temp_coord1 = new float[3];
                    unsigned long n_connect = 0;
                    unsigned long n_coords = 0;
                    for(unsigned long i = 0; i <  cursor; i ++)
                    {
                        for (int ii = 0 ; ii < nspace;ii++)
                            temp_coord1[ii] = coords[n_coords++];
                      
                        pcoords->SetTuple3(i,(double)temp_coord1[0],(double)temp_coord1[1],(double)temp_coord1[2]);
                    }
                    
                    out_ds1->Allocate(totalnumberofcellsinthisobj*cellpoints);// this step is important
                    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
                    points->SetData(pcoords);
                    
                    vtkSmartPointer<vtkFloatArray> pdata = vtkSmartPointer<vtkFloatArray>::New();
                    vtkSmartPointer<vtkFloatArray> pdata1 = vtkSmartPointer<vtkFloatArray>::New();
                    vtkSmartPointer<vtkFloatArray> pdata2 = vtkSmartPointer<vtkFloatArray>::New();
                    vtkSmartPointer<vtkFloatArray> pdata3 = vtkSmartPointer<vtkFloatArray>::New();
                    vtkSmartPointer<vtkFloatArray> pdata4 = vtkSmartPointer<vtkFloatArray>::New();
                    
                    pdata->SetNumberOfTuples((vtkIdType)cursor);
                    pdata->SetName("OW");
                    pdata1->SetNumberOfTuples((vtkIdType)cursor);
                    pdata1->SetName("Temperature");
                    pdata2->SetNumberOfTuples((vtkIdType)cursor);
                    pdata2->SetName("Salt");
                    pdata3->SetNumberOfTuples((vtkIdType)cursor);
                    pdata3->SetName("U");
                    pdata4->SetNumberOfTuples((vtkIdType)cursor);
                    pdata4->SetName("V");
                    
                    double *temp_data = new double;
                    double *temp_data1 = new double;
                    double *temp_data2 = new double;
                    double *temp_data3 = new double;
                    double *temp_data4 = new double;
                    
                    for (unsigned long i = 0; i < cursor; i++)
                    {
                        *temp_data = nodevals[i];
                        *temp_data1 = nodevals1[i];
                        *temp_data2 = nodevals2[i];
                        *temp_data3 = nodevals3[i];
                        *temp_data4 = nodevals4[i];
                        pdata->SetTuple1(i,(double)*temp_data);
                        pdata1->SetTuple1(i,(double)*temp_data1);
                        pdata2->SetTuple1(i,(double)*temp_data2);
                        pdata3->SetTuple1(i,(double)*temp_data3);
                        pdata4->SetTuple1(i,(double)*temp_data4);
                    }                    
                    
                    vtkIdType *temp_connect = new vtkIdType[cellpoints];     //int n_connect = 0;
                    
                    for(unsigned long i = 0; i < totalnumberofcellsinthisobj;i++)
                    {
                        for (int j = 0; j < cellpoints;j++)
                            temp_connect[j] = connects[n_connect++];

                        out_ds1->InsertNextCell(VTK_HEXAHEDRON, (vtkIdType) cellpoints, temp_connect);       //this is for Enrique's data which is saved as structure grid
                        //out_ds1->InsertNextCell(VTK_VOXEL, (vtkIdType) cellpoints, temp_connect);
                    }
                   
                    out_ds1->SetPoints(points);
                    //out_ds1->GetPointData()->SetScalars(pdata);

                    out_ds1->GetPointData()->AddArray(pdata);
                    out_ds1->GetPointData()->AddArray(pdata1); 
                    out_ds1->GetPointData()->AddArray(pdata2); 
                    out_ds1->GetPointData()->AddArray(pdata3); 
                    out_ds1->GetPointData()->AddArray(pdata4); 

                    
                    vtkSmartPointer<vtkContourFilter> isosurface = vtkSmartPointer<vtkContourFilter>::New();
                    isosurface->SetInputData(out_ds1);
                    //isosurface->SetNumberOfContours(1);


                    isosurface->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS , "OW");                    
                    isosurface->SetValue(0, thresh1);
                    isosurface->ComputeGradientsOn();
                    isosurface->ComputeNormalsOn();	

                    isosurface->Update();
                    vtkSmartPointer<vtkPolyData> outpoly2 = vtkSmartPointer<vtkPolyData>::New();
                    outpoly2 = isosurface->GetOutput();
                    //outpoly2->Update();

    
                    delete temp_data;
                    delete temp_data1;
                    delete temp_data2;
                    delete temp_data3;
                    delete temp_data4;
                    delete[] connects;
                    delete[] coords;
                    delete[] nodevals;
                    delete[] nodevals1;
                    delete[] nodevals2;
                    delete[] nodevals3;
                    delete[] nodevals4;
                    delete[] temp_coord1;
                    delete[] temp_connect;
                    
                    
                    
                    double p = 0;
                    ofstream pfile;
                    if (firstobj ==0)
                       pfile.open(mypolyfile.c_str(),ofstream::out | ofstream::trunc);
             
                    else
                       pfile.open(mypolyfile.c_str(),ofstream::out | ofstream::app);
                 
                    if (!pfile.is_open())
                    {
                       cout<<"ExtractSurf:cannot open "<<mypolyfile<<endl;
                       return 0;
                    }
                    
                    
                    vtkSmartPointer<vtkIdList> listconnect = vtkSmartPointer<vtkIdList>::New();

                    vtkSmartPointer<vtkFloatArray> pointsvalue_OW = vtkFloatArray::SafeDownCast(outpoly2->GetPointData()->GetArray("OW"));
                    vtkSmartPointer<vtkFloatArray> pointsvalue_Temp = vtkFloatArray::SafeDownCast(outpoly2->GetPointData()->GetArray("Temperature"));
                    vtkSmartPointer<vtkFloatArray> pointsvalue_Salt = vtkFloatArray::SafeDownCast(outpoly2->GetPointData()->GetArray("Salt"));
                    vtkSmartPointer<vtkFloatArray> pointsvalue_U = vtkFloatArray::SafeDownCast(outpoly2->GetPointData()->GetArray("U"));
                    vtkSmartPointer<vtkFloatArray> pointsvalue_V = vtkFloatArray::SafeDownCast(outpoly2->GetPointData()->GetArray("V"));

                    string frame_date = computedate(currentTime+1);
               
                    pfile<<"******day******"<<endl;
                    pfile<<frame_date<<endl;
                    pfile<<155<<" "<<155<<" "<<155<<endl;
                    pfile<<setprecision(5);
                    pfile<<spin_orientation<<endl;
                    pfile<< volume << " " << Cx << " " << Cy << " " << Cz << " " << LowerLeft_Corner_x << " " << LowerLeft_Corner_y << " " << LowerLeft_Corner_z << " " << UpperRight_Corner_x << " " << UpperRight_Corner_y << " " << UpperRight_Corner_z << endl;
                    pfile<< outpoly2->GetNumberOfPoints() << " " << MinTemp << " " << MaxTemp << " " << MinSalt << " " << MaxSalt << " " << sqrt(MinSpeed) << " " << sqrt(MaxSpeed) << endl;
                    pfile.setf(ios::fixed,ios::floatfield);
                    pfile<<setprecision(3);
                    double *temp_coord_array = new double[3];

                    for(long k3=0; k3<outpoly2->GetNumberOfPoints(); k3++)
                    {
                        temp_coord_array = outpoly2->GetPoint(k3);
                        for (int k4 = 0; k4 <nspace;k4++)
                            pfile<<temp_coord_array[k4] << " ";
                       
                        //pfile << pointsvalue_Temp->GetValue(k3) << " " << pointsvalue_Salt->GetValue(k3) << " " << pointsvalue_U->GetValue(k3) << " " << pointsvalue_V->GetValue(k3) << " ";
                        
                        pfile << pointsvalue_Temp->GetValue(k3) << " " << pointsvalue_Salt->GetValue(k3) << " " << sqrt( pointsvalue_U->GetValue(k3)*pointsvalue_U->GetValue(k3) + pointsvalue_V->GetValue(k3)*pointsvalue_V->GetValue(k3) ) << " ";
    
                        
                        pfile<<endl;
                    }
                    //delete[] temp_coord_array;
                    
                    pfile<<outpoly2->GetNumberOfCells()<<endl;
                    for(long i=0;i<outpoly2->GetNumberOfCells();i++)
                    {
                        outpoly2->GetCellPoints(i,listconnect);
                        pfile<<"3 ";
                        for(long j = 0; j < listconnect->GetNumberOfIds();j++)
                        {
                            unsigned long temp_id = (unsigned long)listconnect->GetId((const vtkIdType) j) ;
                            pfile<< temp_id << " ";
                        }
                        pfile << endl;
                    }
                    pfile<<0<<endl<<endl;
                    //#cout << "finished writing poly file\n";
                    pfile.close();
 
                }
            }
            
        }
        
        cout<<"Eddy Extraction is completed.. Total number of extracted eddies in this time step is: [" <<firstobj+1<<"]" <<endl;
        
        if (firstobj+1 == 0)
            cout<<"---- WARNING:  No Object is found with the given parameters!  --- "<<endl;
        

        
        in_ds->Delete();
        in_ds = NULL;
        
        char mypolyfile1[256];
        //sprintf(mypolyfile1, mypolyfile.c_str());
        
        sprintf(mypolyfile1, "%s", mypolyfile.c_str());
 
        
        int finished = 0;
        
        cout<<"The Directory of Output Files: "<<base_GeneratedTrackFileName << endl;
        cout<< "Start Eddy Tracking !! "<<endl;
        
        int a = BeginFeatureTrack(base_GeneratedTrackFileName, currenttimevalue, mypolyfile1, finished, listfile, InitialtimeStep ,TimeIncrement, currentTime, TimePrecision);
        
        cout<< "Finish Eddy Tracking !! "<<endl;
        
        
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
 
        cout<< endl;
        cout << "**************** Time elapsed: " << elapsed_secs << " seconds for computing the "<<currentTime+1<<"th time step!"<< " ****************" << endl;
        cout<< endl;

        delete[] ProcessedCell;
        
        remove(zetafile.c_str());
        
        
    } // end of currenttime loop


    
}//end of main




