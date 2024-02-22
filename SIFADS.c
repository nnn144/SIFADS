/** ==============================================================================
 ** ======  SIFADS (Spline Initialized FADS) algorithm for dynamic SPECT =========
 ** ==============================================================================
 **
 **     **********************   Version 1.2.0         **********************
 **
 **     **********************   Copyright (c)  2016   **********************
 **
 ** --------------------------------------------------------------------------------
 ** Author:       Mahmoud Abdalah
 **
 ** Institute:    Florida Institute of Technology & Lawrence National Berkeley Lab
 ** --------------------------------------------------------------------------------
 **
 **
 ** This code is the implementation of SIFADS algorithm to estimate time activity
 ** curves directly from projection data (dynamic reconstruction). A detailed description
 ** of SIFADS algorithm can be found in the following paper:
 **
 **            Abdalah M, Boutchko R, Mitra D and Gullberg G.T. "Reconstruction of 4-D
 **            Dynamic SPECT Images From Inconsistent Projections Using a Spline Initialized
 **            FADS Algorithm (SIFADS)". IEEE Trans Med Imaging,34(1): 216-228, 2015.
 **
 ** Also, the code preforms static SPECT reconstruction, back projection, and forward projection.
 **
 **
 **
 **
 ** To compile:
 **               gcc -o SIFADS SIFADS.c
 ** To run:
 **                 ./SIFADS {argument(1), argument(2), argument(3), ..., argument(n)}
 **
 ** Run the code with empty arguments for help.
 **
 ** Examples of input arguments and output files:
 **
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **                                                 Static forward projection
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "F";                   //Operation
 **    argv[2] = "SystemMatrix.dat";    //Name of system matrix file                                ForwardProjected.Static.Sinogram.sino
 **    argv[3] = "Volume.vol";          //Name of volume file
 **    argv[4] = "Sinogram.sino";       //Name of sinogram file
 **
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **                                                   Static back projection
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "B";                   //Operation
 **    argv[2] = "SystemMatrix.dat";    //Name of system matrix file                                 backprojected.Static.Volume.vol
 **    argv[3] = "Volume.vol";          //Name of volume file
 **    argv[4] = "Sinogram.sino";       //Name of sinogram file
 **
 ** ---------------------------------------------------------------------------------------------------------------------------------------
 **                                                  Static reconstruction
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "R";                   //Operation
 **    argv[2] = "SystemMatrix.dat";    //Name of system matrix file                                  Reconstructed.Static.Volume.vol
 **    argv[3] = "Volume.vol";          //Name of volume file
 **    argv[4] = "Sinogram.sino";       //Name of sinogram file
 **    argv[5] = “64”;                  //Size of volume in x direction
 **    argv[6] = “64”;                  //Size of volume in y direction
 **    argv[7] = “64”;                  //Size of volume in z direction
 **
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **                                               Dynamic forward projection
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "TF";                   //Operation
 **    argv[2] = "SystemMatrix.dat";     //Name of system matrix file
 **    argv[3] = "Coefficients.coff";    //Name of coefficients file                             ForwardProjected.Dynamic.Sinogram.sino
 **    argv[4] = "TACs.txt";             //Name of TACs file
 **    argv[5] = "Sinogram.sino";        //Name of sinogram file
 **    argv[6] = "N_ROTATIONS";
 **    argv[7] = "4096";
 **    argv[8] = “2”;
 **
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **                                                 Dynamic back projection
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "TB";                   //Operation
 **    argv[2] = "SystemMatrix.dat";     //Name of system matrix file
 **    argv[3] = "Coefficients.coff";    //Name of coefficients file                              BackProjected.Dynamic.Coefficients.coff
 **    argv[4] = "TACs.txt";             //Name of TACs file
 **    argv[5] = "Sinogram.sino";        //Name of sinogram file
 **    argv[6] = “1”;                    //No of angular rotations in the sinogram
 **    argv[7] = "4096";                 //No of bin on the head/detector
 **    argv[8] = “2”;                    //No of heads/detectors
 **
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **                                             Dynamic reconstruction (SIFADS)
 **  ---------------------------------------------------------------------------------------------------------------------------------------
 **  Inputs:                                                                                     Outputs:
 **  ------                                                                                      --------
 **    argv[1] = "TR";                   //Operation
 **    argv[2] = "SystemMatrix.dat";     //Name of system matrix file                               Final_Factors.csv
 **    argv[3] = "Coefficients.coff";    //Name of coefficients file                                Final_Coefficients.coff
 **    argv[4] = "Bsplines.txt";         //Name of splines file                                     Final_TACs.csv
 **    argv[5] = "Sinogram.sino";        //Name of sinogram file
 **    argv[6] = “1”;                    //No of angular rotations in the sinogram
 **    argv[7] = “4096”;                 //No of bin on the head/detector
 **    argv[8] = “2”;                    //No of heads/detectors
 **    argv[9] = "Segmented_volume.vol"; //Name of segmented static volume file
 **    argv[10] = “3”;                   //No of segments in the segmented static volume
 **    argv[11] = “64”;                  //Size of volume in x direction
 **    argv[12] = “64”;                  //Size of volume in y direction
 **    argv[13] = “64”;                  //Size of volume in z direction
 ** =======================================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//==================================================================================
// ========================= parameters that may be adjusted      ==================
//==================================================================================
// input sinogram can be either 16-bit unsigned integer (unsigned short) or float
// uncomment this, if input sinogram is float,comment if 16-bit unsigned integer
//#define FLOAT_INPUT_SINOGRAM

// Volume is best output as 4-bit float
#define FLOAT_INPUT_VOLUME

// No of iterations
// Typically MLEM needs around 20 iterations
#define N_ITERATIONS_MLEM 10       // MLEM max iterations

// Provide the needed values for the regularization parameters.
#define DEFAULT_SMOOTHNESS_LAMBDA 0.001      // Total variation lambda
#define DEFAULT_MIX_LAMBDA 100               // coefficient mix lambda
#define DEFAULT_TB_SMOOTHNESS_LAMBDA 0.9      // Curve smoothing lambda

// Uncomment this, if need to use static mask (segmented static image) in reconstruction
// comment if not using the static mask in reconstruction
#define STATIC_MASKING

// Do not change these parameters
// This defines how small zero is, and what is one divided by zero
#define ONE_OVER_EPSILON 3.0e+13
#define EPSILON 3.0e-13

size_t size_float= sizeof(float);
size_t size_int= sizeof(int);

// a structure to store all passed parameters between functions.
// this is just for code simplisty.
struct PARAMETERS
{
    float *sinogram;           // Sinogram array
    char *sysmat_fname;        // name of system matrix
    float *StaticMask;         // mask created out of segmented static reconstruction
    
    int N_rotations;           // Number of angular rotations included in the input sinogram
    int NumBins;  // size of a projection/number of pins on head
    int basisNum;              // number of time basis/factors/TACs
    int basisLen;              // length of time basis/factors/TACs
    int Nvolume;               // size of imaged volume
    int Ncoefficients;         // size of coefficients = Nvolume * basisNum
    int Nsinogram;             // size of sinogram
    int No_Heads;              // number of detectors/heads use for acquiring the data
    int Nx;                    // dimenrsion of imaged volume in x direction
    int Ny;                    // dimenrsion of imaged volume in y direction
    int Nz;                    // dimenrsion of imaged volume in z direction
    
    float TB_SMOOTHNESS_LAMBDA; // Regularization parameter of curves Smoothness
    float SMOOTHNESS_LAMBDA;    // Regularization parameter of coefficient smoothness (total variation)
    float MIX_LAMBDA;           // Regularization parameter of coefficient mix
};

// returns the voxel position in 1d array by given x,y,z coords and,
// c= total number of cols in x direction and,
// r= total number of rows in y direction
#define voxel_pos(r,c, x,y,z) (x+y*r+z*r*c)
// calculates the gamma value for dynamic lambda estimation
#define GAMMA(a,b,c) (c*(pow((a/(0.05*b)),0.25)))

// --- functions prototypes:
int main(int argc, char **argv);
//=========================================
void Do_static_BackProjection(char **argv);
void Do_static_ForProjection(char **argv);
void Do_static_Recon(char **argv);

void Do_dynamic_BackProjection(char **argv);
void Do_dynamic_ForProjection(char **argv);
void Do_dynamic_Recon(char **argv);
//=========================================
void ST_MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, float *volume,int Nx,
                   int Ny, int Nz,float SMOOTHNESS_LAMBDA, char *fname );
void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,char *fname);
void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,char *fname);
//=============================================================================================
void Spline_MLEM_recon(float *coeffVolume,float *basis, int N_iter,struct PARAMETERS *func_parameters);
void FADS_MLEM_recon(float *coeffVolume,float *factors,struct PARAMETERS *func_parameters);
void TimeBasisBP_c(int Nsinogram,int No_Heads, int N_rotations, int Ncoefficients, int
                   NumBins,int basisNum, int basisLen, float *sinogram, float
                   *coeffVolume, float *denominator_c,float *basis, char *sysmat_fname);
void TimeBasisBP_f(int Nsinogram,int No_Heads, int N_rotations, int Ncoefficients, int
                   NumBins, int basisNum, int basisLen, float *sinogram,
                   float *coeffVolume, float *denominator_f, float *basis, char *sysmat_fname);
void TimeBasisFP(int Nsinogram,int No_Heads, int N_rotations, int Ncoefficients,
                 int NumBins,int basisNum, int basisLen, float *sinogram, float
                 *coeffVolume,float *basis, char *sysmat_fname);

//===============================================================================================
float NN_3D(float *NN_grad, float *volume,float *mask, int Nx, int Ny, int Nz, int omega);
float Coeff_Mix(float *MixGrad, float *coefficients,int Nvolume, int basisNum,float *mask);
float TB_Smooth(float *TBSmooth_grad,float *basis, int basisLen, int basisNum);
//===============================================================================================
void Convert2Mask(float *mask, float *volume,int Nvolume,float Threshold);
void CreateMask(float *Dynamic_Mask, float *Static_Mask, float *coefficients,
                int Nvolume,int basisNum,struct PARAMETERS *func_parameters);
void convolve_coeffs(float *Ave_NN, float *coefficients,float *mask, int Nvolume,struct
                     PARAMETERS *func_parameters);
float otsu(float *volume, int size);
//===============================================================================================
void forwardDFT_1D(const float *data, const int N, float *real, float *imaginary);
void inverseDFT_1D(const float *real, const float *imaginary, const int N, float *data);
void Smooth_curve( float *curve, int basisLen);
//===============================================================================================

void Normailze_coefficients(float *coefficients, float *basis, int basisNum, int tbLen,
                            int Nvolume);
void Compute_TACS_Volume(float *Tvolumes,float *volume, float *basis, int tbNum, int tbLen,
                         int Nvolume);
void Compute_Averaged_TACS(float *Tvolumes,float *segmentedVolume, float *basis, int tbNum,
                           int tbLen,int Nvolume);
void Read_volume(float *Data, int read_size,int Data_Type,  char *fname);
void Save_Volume(float *volume, int Save_size,char Save_Name[128],char Save_Extension[128]);
void Read_Curves(float *basis, int tbNum,int tbLen,  char *basis_fname);
void Print_Curves(float *basis, int tbNum, int tbLen, char title[128]);
void Save_Curves(float *basis, int tbNum, int tbLen, char Save_Name[128]);
double getCPUTime();
//===============================================================================================
// ----------------------------------- END HEADER, BEGIN FUNCTIONS ------------------------------

// --- main code.
//     command line input arguments are explained in error message printed when a.out is started proj_weights
int main(int argc, char **argv) {
    int tbFlag = 0;
    char Action, *opArg;
    
    if  (argc == 1)
    {
        fprintf(stderr,"\n-----------------------------------------------------\n");
        fprintf(stderr,"STATIC USAGE:\n");
        
        fprintf(
                stderr,
                "Static reconstruction arguments: \n\nOperation Sysmat.file Volume.file Sinogram.file  Volume_x Volume_y Volume_z \n\n");
        fprintf(
                stderr,
                "Where \"Operation\" is R for reconstruction, F for forward projection, or B for back projection. ");
        fprintf(
                stderr,"Volume_x, Volume_y, and Volume_z are the volume dimensions (Dimensions are needed in case of TV is used in the static reconstruction)\n\n");
        
        fprintf(stderr,"\n-----------------------------------------------------\n");
        fprintf(stderr,"DYNAMIC USAGE:\n");
        
        fprintf(
                stderr,
                "Dynamic reconstruction arguments:\n\nOperation Sysmat.file Coefficients.file Timebasis.file Sinogram.file N_rotations NumBins N_heads Segmented.Static.Volume NumnerOfSegments Volume_x Volume_y Volume_z\n\n");
        fprintf(
                stderr,
                "Where \"operation\" is TR for reconstruction, TF for forward projection, TB for forward projection. Splines.file is text file contains the the set of b-splines used for dynamic reconstruction.");
        fprintf(
                stderr,
                "N_rotations is the number of of angular rotations included in the input sinogram. NumBins is the number of bins on the detector. NumberOfHeads can be 1 or 2. NumnerOfSegments is the number of segments in the Segmented.Static.Volume. Volume_x, Volume_y, and Volume_z are the static volume dimensions (Dimensions are always needed in dynamic reconstruction). \n");
        exit(1);
    }
    
    
    // find the entered action
    Action = ' ';
    opArg =strpbrk(argv[1], "rRfFbB");
    
    if (opArg == NULL)
    {
        fprintf(
                stderr,
                "Action setting should be: \"R\" or \"r\"  for reconstruction, \"F\" or \"f\" for forwardprojection,\n");
        fprintf(
                stderr,
                "or \"B\" or \"b\" for backprojection, unrecongnized action \"%s\"\n",
                argv[1]);
        exit(1);
    }
    else
    {
        switch (*opArg) {
            case 'r':
            case 'R':
                Action = 'r';
                break;
                
            case 'f':
            case 'F':
                Action = 'f';
                break;
                
            case 'b':
            case 'B':
                Action = 'b';
                break;
        }
    }
    
    // Check if the action is for dynamic reconstruction
    opArg = strpbrk(argv[1], "Tt");
    if (opArg != NULL) {
        tbFlag = 1;
    }
    
    //  --- Main action ---
    switch (Action) {
        case 'r':
            if (tbFlag)  // Dynamic reconstruction
            {
                if (argc != 14)
                {
                    fprintf(
                            stderr,
                            "Dynamic reconstruction arguments:\n\nOperation Sysmat.file Coefficients.file Splines.file Sinogram.file N_rotations NumBins N_heads Segmented.Static.Volume NumnerOfSegments Volume_x Volume_y Volume_z\n\n");
                    fprintf(
                            stderr,
                            "Where \"operation\" is TR for reconstruction, TF for forward projection, TB for forward projection. Splines.file is text file contains the the set of b-splines used for dynamic reconstruction.");
                    fprintf(
                            stderr,
                            "N_rotations is the number of of angular rotations included in the input sinogram. NumBins is the number of bins on the detector. NumberOfHeads can be 1 or 2. NumnerOfSegments is the number of segments in the Segmented.Static.Volume. Volume_x, Volume_y, and Volume_z are the static volume dimensions (Dimensions are always needed in dynamic reconstruction). \n");
                    exit(1);
                }
                
                Do_dynamic_Recon(argv);
            }
            else // Static reconstruction
            {
                if  (argc != 8)
                {
                    fprintf(
                            stderr,
                            "Static reconstruction arguments: \n\nOperation Sysmat.file Volume.file Sinogram.file  Volume_x Volume_y Volume_z \n\n");
                    fprintf(
                            stderr,
                            "Where \"Operation\" is R for reconstruction, F for forward projection, or B for back projection. ");
                    fprintf(
                            stderr,"Volume_x, Volume_y, and Volume_z are the volume dimensions (Dimensions are needed in case of TV is used in the static reconstruction)\n\n");
                    
                    exit(1);
                }
                
                Do_static_Recon(argv);
            }
            break;
        case 'f':
            if (tbFlag)  // Dynamic forward projection
            {
                if  (argc != 9)
                {
                    fprintf(
                            stderr,
                            "\nDynamic forward projection:\n TF Sysmat.file Coefficients.file TACs.file Sinogram.file N_rotations NumBins N_heads. \n");
                    fprintf(
                            stderr,
                            "Where: N_rotations is the number of of angular rotations included in the input sinogram. NumBins is the number of bins on the detector. NumberOfHeads can be 1 or 2. TACs.file is a text file contains the the set of TACs functions used to describe the temporal activity in each tissue.\n");
                    exit(1);
                }
                
                Do_dynamic_ForProjection(argv);
            }
            else // Static forward projection
            {
                if  (argc != 5)
                {
                    fprintf(
                            stderr,
                            "\nStatic forward projection:\n F Sysmat.file Volume.file Sinogram.file \n");
                    fprintf(
                            stderr,
                            "Where \"operation\" is F for forward.\n");
                    
                    exit(1);
                }
                
                Do_static_ForProjection(argv);
            }
            break;
        case 'b':
            if (tbFlag)  // Dynamic back projection
            {
                if  (argc != 9)
                {
                    fprintf(
                            stderr,
                            "\nDynamic back projection:\n TB Sysmat.file Coefficients.file TACs.file Sinogram.file N_rotations NumBins N_heads. \n");
                    fprintf(
                            stderr,
                            "Where: N_rotations is the number of of angular rotations included in the input sinogram. NumBins is the number of bins on the detector. NumberOfHeads can be 1 or 2. TACs.file is a text file contains the the set of TACs functions used to describe the temporal activity in each tissue.\n");
                    exit(1);
                }
                
                Do_dynamic_BackProjection(argv);
            }
            else // Static back projection
            {
                if  (argc != 5)
                {
                    fprintf(
                            stderr,
                            "\nStatic forward projection:\n B Sysmat.file Volume.file Sinogram.file \n");
                    fprintf(
                            stderr,
                            "Where \"operation\" is B for back projection.\n");
                    
                    exit(1);
                }
                
                Do_static_BackProjection(argv);
            }
            break;
    }
    fprintf(stderr,"\n Done! \n");
}

/*
 * Performs static forward projection
 */
void Do_static_ForProjection(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram;
    float *sinogram, *volume= NULL;
    char sino_fname[128], volume_fname[128], sysmat_fname[128], buf[128];
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(volume_fname, argv[3]);
    strcpy(sino_fname, argv[4]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // allocate sinogram and volume
    sinogram = (float *) malloc(size_float * Nsinogram);
    volume = (float *) malloc(size_float * Nvolume);
    
    //-------------- read volume
    // 1 = 8-bit singram, 2 = 32-bit sinogram (float)
#ifndef FLOAT_INPUT_VOLUME
    Read_volume(volume,  Nvolume,1,  image_fname);
#else
    Read_volume(volume,  Nvolume,2,  volume_fname);
#endif
    
    fprintf( stderr,"Forward project volume of length %d into sinogram of length %d\n",Nvolume, Nsinogram);
    RayDrFP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
    
    // Save final volume
    strcpy(buf, "ForwardProjected.Static.");
    strcat(buf,sino_fname);
    
    Save_Volume(sinogram,Nsinogram,buf,"sino");
    
    free(sinogram);
    free(volume);
    return;
}

/*
 * Performs static back projection
 */
void Do_static_BackProjection(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram;
    float *sinogram, *volume= NULL;
    char sino_fname[128], volume_fname[128], sysmat_fname[128], buf[128];
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(volume_fname, argv[3]);
    strcpy(sino_fname, argv[4]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // allocate sinogram and volume
    sinogram = (float *) malloc(size_float * Nsinogram);
    volume = (float *) malloc(size_float * Nvolume);
    
    
    //-------------- read sinogram
    // 1 = 8-bit singram, 2 = 32-bit sinogram (float)
#ifndef FLOAT_INPUT_SINOGRAM
    
    Read_volume(sinogram,  Nsinogram,1,  sino_fname);
#else
    Read_volume(sinogram,  Nsinogram,2,  sino_fname);
#endif
    
    fprintf(stderr,"Backproject volume of length %d from sinogram of length %d\n",Nvolume, Nsinogram);
    RayDrBP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
    
    // Save final volume
    strcpy(buf, "backprojected.Static.");
    strcat(buf,volume_fname);
    
    Save_Volume(volume,Nvolume,buf,"Vol");
    
    free(sinogram);
    free(volume);
    return;
}

/*
 * Performs static reconstruction
 */
void Do_static_Recon(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram, Nx,Ny,Nz;
    float *sinogram, *volume= NULL;
    char sino_fname[128], volume_fname[128], sysmat_fname[128],buf[128];
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(volume_fname, argv[3]);
    strcpy(sino_fname, argv[4]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // copy the volume dimensions in case of static reconstruction with TV or dynamic reconstruction
    Nx = atoi(argv[5]);
    Ny =atoi(argv[6]) ;
    Nz = atoi(argv[7]);
    
    // check if the size of volume is consistent with the input dimensions
    if(Nvolume != (Nx*Ny*Nz))
    {
        fprintf(
                stderr,
                "The input volume dimensions %dx%dx%d = %d are not consistent with the volume size in the system matrix %d \n",Nx,Ny,Nz,Nx*Ny*Nz,Nvolume);
        exit(1);
    }
    
    // allocate sinogram and volume
    sinogram = (float *) malloc(size_float * Nsinogram);
    volume = (float *) malloc(size_float * Nvolume);
    
    // check if the size of volume is consistent with the input dimensions
    if(Nvolume != (Nx*Ny*Nz))
    {
        fprintf(
                stderr,
                "The input volume dimensions are not consistent with the volume size in the system matrix \n");
        exit(1);
    }
    
    //-------------- read sinogram
#ifndef FLOAT_INPUT_SINOGRAM
    
    // 1 = 8-bit singram, 2 = 32-bit sinogram (float)
    Read_volume(sinogram,  Nsinogram,1,  sino_fname);
#else
    Read_volume(sinogram,  Nsinogram,2,  sino_fname);
#endif
    
    fprintf( stderr,"Static MLEM reconstruction volume of length %d from sinogram of length %d\n",Nvolume, Nsinogram);
    ST_MLEM_recon(Nsinogram, sinogram, Nvolume, volume, Nx, Ny, Nz,DEFAULT_SMOOTHNESS_LAMBDA,sysmat_fname);
    
    // Save final volume
    strcpy(buf, "Reconstructed.Static.");
    strcat(buf,volume_fname);
    
    Save_Volume(volume,Nvolume,buf,"Vol");
    
    free(sinogram);
    free(volume);
    return;
}

/*
 * Performs dynamic forward projection
 */
void Do_dynamic_BackProjection(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram, Ncoefficients,No_Heads;
    int tbNum, tbLen,NumBins, N_rotations;
    float *sinogram, *coefficients, *TACs, *tmp_coefficients;
    char sino_fname[128], coeffs_fname[128], sysmat_fname[128],
    TACs_fname[128],buf[128];
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(coeffs_fname, argv[3]);
    strcpy(TACs_fname, argv[4]);
    strcpy(sino_fname, argv[5]);
    
    // Copy the numebr of rotations
    N_rotations=atoi(argv[6]);
    
    // size of the head/projection
    NumBins = atoi(argv[7]);
    
    // Copy the numebr of heads
    No_Heads=atoi(argv[8]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // first read the number and length of time basis
    if ((fid = fopen(TACs_fname, "r")) == NULL) {
        fprintf(stderr, "Could not open basis file %s for reading.\n",
                TACs_fname);
        exit(1);
    }
    fscanf(fid, "%d\n", &tbNum);
    fscanf(fid, "%d\n", &tbLen);
    fclose(fid);
    
    // read the contents of the time basis file into an array
    Read_Curves(TACs, tbNum, tbLen, TACs_fname);
    
    // check if the entered sinogram is consistent with sinogram in the system matrix
    if ((NumBins <= 0) || (NumBins > Nsinogram) || (NumBins*tbLen*No_Heads != Nsinogram)) {
        fprintf(stderr,"Number of bins %d must be an integer smaller than the size of the original sinogram %d.\n",
                NumBins, Nsinogram);
        exit(1);
    }
    fprintf(stderr, "Using %d basis functions of length %d.\n", tbNum,
            tbLen);
    
    // define the size of coefficients
    Ncoefficients = Nvolume * tbNum;
    
    // allocate sinogram and coefficients arrays
    sinogram = (float *) malloc(size_float * Nsinogram);
    coefficients = (float *) malloc(size_float * Ncoefficients);
    
    //-------------- read sinogram
    // 1 = 8-bit singram, 2 = 32-bit sinogram (float)
#ifndef FLOAT_INPUT_SINOGRAM
    Read_volume(sinogram,  Nsinogram,1,  sino_fname);
#else
    Read_volume(sinogram,  Nsinogram,2,  sino_fname);
#endif
    
    fprintf(stderr,"Backproject coefficients of length %d from sinogram of length %d\n",Ncoefficients, Nsinogram);
    
    tmp_coefficients = (float *) malloc(size_float * Ncoefficients);
    
    TimeBasisBP_c(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, tbNum, tbLen,
                  sinogram, coefficients,tmp_coefficients, TACs, sysmat_fname);
    
    // Save final volume
    strcpy(buf, "Dynamic_backprojection.Dynamic.");
    strcat(buf,coeffs_fname);
    
    Save_Volume(coefficients,Ncoefficients,buf,"coff");
    
    free(tmp_coefficients);
    free(coefficients);
    free(sinogram);
    free(TACs);
    return;
}

/*
 * Performs dynamic back projection
 */
void Do_dynamic_ForProjection(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram, N_rotations, Ncoefficients,No_Heads;
    int tbNum, tbLen,NumBins;
    float *sinogram, *coefficients,*TACs;
    char sino_fname[128], coeffs_fname[128], sysmat_fname[128],
    TACs_fname[128],buf[128];
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(coeffs_fname, argv[3]);
    strcpy(TACs_fname, argv[4]);
    strcpy(sino_fname, argv[5]);
    
    // Copy the numebr of rotations
    N_rotations=atoi(argv[6]);
    
    // size of the head/projection
    NumBins = atoi(argv[7]);
    
    // Copy the numebr of heads
    No_Heads=atoi(argv[8]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // first read the number and length of time basis
    if ((fid = fopen(TACs_fname, "r")) == NULL) {
        fprintf(stderr, "Could not open basis file %s for reading.\n",
                TACs_fname);
        exit(1);
    }
    fscanf(fid, "%d\n", &tbNum);
    fscanf(fid, "%d\n", &tbLen);
    fclose(fid);
    
    // read the contents of the time basis file into an array
    Read_Curves(TACs, tbNum, tbLen, TACs_fname);
    
    // check if the entered sinogram is consistent with sinogram in the system matrix
    if ((NumBins <= 0) || (NumBins > Nsinogram) || (NumBins*tbLen*No_Heads != Nsinogram)) {
        fprintf(stderr,"Number of bins %d must be an integer smaller than the size of the original sinogram %d.\n",
                NumBins, Nsinogram);
        exit(1);
    }
    fprintf(stderr, "Using %d basis functions of length %d.\n", tbNum,
            tbLen);
    
    // define the size of coefficients
    Ncoefficients = Nvolume * tbNum;
    
    // allocate sinogram and coefficients arrays
    sinogram = (float *) malloc(size_float * Nsinogram);
    coefficients = (float *) malloc(size_float * Ncoefficients);
    
    //-------------- read coefficients
    // 1 = 8-bit coefficients, 2 = 32-bit coefficients (float)
#ifndef FLOAT_INPUT_VOLUME
    Read_volume(coefficients,  Ncoefficients,1,  coeffs_fname);
#else
    Read_volume(coefficients,  Ncoefficients,2,  coeffs_fname);
#endif
    
    fprintf( stderr,"Forward project coefficiets of length %d into sinogram of length %d\n",Ncoefficients, Nsinogram);
    
    TimeBasisFP(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, tbNum, tbLen,
                sinogram, coefficients, TACs, sysmat_fname);
    
    
    // Save final volume
    strcpy(buf, "Dynamic_forwardprojection.Dynamic.");
    strcat(buf,sino_fname);
    
    Save_Volume(coefficients,Ncoefficients,sino_fname,"sino");
    
    free(coefficients);
    free(sinogram);
    free(TACs);
    return;
}

/*
 * Performs dynamic reconstruction (SIFADS)
 */
void Do_dynamic_Recon(char **argv)
{
    FILE *fid;
    int Nvolume, Nsinogram, Nx,Ny,Nz,Ncoefficients,No_Heads,N_rotations;
    int tbNum, tbLen, NumBins, numbSegments;
    float *Splines,*factors,*Tvolumes;
    float *sinogram, *coefficients;
    float *StaticMask;
    struct PARAMETERS func_parameters;
    char sino_fname[128], coeffs_fname[128], sysmat_fname[128],Segmented_Volume_fname[128],
    Spline_fname[128],buf[128];
    
    
    // Copy files' names
    strcpy(sysmat_fname, argv[2]);
    strcpy(coeffs_fname, argv[3]);
    strcpy(Spline_fname, argv[4]);
    strcpy(sino_fname, argv[5]);
    strcpy(Segmented_Volume_fname, argv[9]);
    
    Nx = atoi(argv[11]);
    Ny =atoi(argv[12]) ;
    Nz = atoi(argv[13]);
    
    // number of rotations
    N_rotations=atoi(argv[6]);
    
    // size of the head/projection
    NumBins = atoi(argv[7]);
    
    // Copy the numebr of heads
    No_Heads=atoi(argv[8]);
    
    // number of segemets provided by us
    numbSegments = atoi(argv[10]);
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // check if the size of volume is consistent with the input dimensions
    if(Nvolume != (Nx*Ny*Nz))
    {
        fprintf(
                stderr,
                "The input volume dimensions %dx%dx%d = %d are not consistent with the volume size in the system matrix %d \n",Nx,Ny,Nz,Nx*Ny*Nz,Nvolume);
        exit(1);
    }
    
    // first read the number and length of splines
    if ((fid = fopen(Spline_fname, "r")) == NULL) {
        fprintf(stderr, "Could not open basis file %s for reading.\n",
                Spline_fname);
        exit(1);
    }
    fscanf(fid, "%d\n", &tbNum);
    fscanf(fid, "%d\n", &tbLen);
    fclose(fid);
    
    // declare the spline array
    Splines = (float *) malloc(size_float * tbLen * tbNum);
    
    // read the contents of the spline file into an array
    Read_Curves(Splines, tbNum, tbLen, Spline_fname);
    
    //multiply the zise of sinogram by the number of rotation used to do the reconstruction
    Nsinogram *= N_rotations;
    
    
    // check if the entered sinogram is consistent with sinogram in the system matrix
    if ((NumBins <= 0) || (NumBins > Nsinogram) || (NumBins*tbLen*No_Heads != Nsinogram)) {
        fprintf(stderr,"Number of bins %d must be an integer smaller than the size of the original sinogram %d.\n",
                NumBins, Nsinogram);
        exit(1);
    }
    fprintf(stderr, "Using %d splines of length %d.\n", tbNum,tbLen);
    
    // define the size of coefficients
    Ncoefficients = Nvolume * tbNum;
    
    // allocate sinogram and coefficients arrays
    sinogram = (float *) malloc(size_float * Nsinogram);
    coefficients = (float *) malloc(size_float * Ncoefficients);
    
    //define and read segments
    StaticMask = (float *) malloc(size_float * numbSegments * Nvolume);
    
    // time-dependent volumes
    Tvolumes = (float *) malloc(size_float * Nvolume*tbLen);
    
    //read segments into mask
    // 1 = 16-bit volume, 2 = 32-bit volume (float)
    Read_volume(StaticMask, (numbSegments) * Nvolume, 2, Segmented_Volume_fname);
    
    //convert segmented volume to a binary mask
    Convert2Mask(StaticMask, StaticMask,(numbSegments) * Nvolume,EPSILON);
    
    // Copy the parameters for dynamic case to the parameters structure
    // function parameters will be passed as stucture of needed parameters.
    func_parameters.sinogram=sinogram;
    func_parameters.sysmat_fname= sysmat_fname;
    func_parameters.Nvolume= Nvolume;
    func_parameters.Nsinogram= Nsinogram;
    func_parameters.N_rotations = N_rotations;
    func_parameters.Ncoefficients= Ncoefficients;
    func_parameters.Nx= Nx;
    func_parameters.Ny= Ny;
    func_parameters.Nz= Nz;
    func_parameters.NumBins= NumBins;
    func_parameters.basisLen= tbLen;
    func_parameters.basisNum= tbNum;
    func_parameters.No_Heads= No_Heads;
    func_parameters.StaticMask= StaticMask;
    func_parameters.SMOOTHNESS_LAMBDA = DEFAULT_SMOOTHNESS_LAMBDA;
    func_parameters.MIX_LAMBDA =DEFAULT_MIX_LAMBDA;
    func_parameters.TB_SMOOTHNESS_LAMBDA= DEFAULT_TB_SMOOTHNESS_LAMBDA;
    
    //-------------- read sinogram
    // 1 = 8-bit singram, 2 = 32-bit sinogram (float)
#ifndef FLOAT_INPUT_SINOGRAM
    Read_volume(sinogram,  Nsinogram,1,  sino_fname);
#else
    Read_volume(sinogram,  Nsinogram,2,  sino_fname);
#endif
    
    //========================================================== B-spline method===============================================
    
    fprintf(stderr, "=============  B-splines Fitting  =========\n");
    // small smoothing in spline method
    func_parameters.SMOOTHNESS_LAMBDA =.001;
    // turn mix regularization function off
    func_parameters.MIX_LAMBDA = 0.0;
    
    fprintf( stderr,"MLEM Reconstruction coefficients of length %d from sinogram of length %d\n",Ncoefficients, Nsinogram);
    
    // reconstruct coefficients by MLEM
    Spline_MLEM_recon(coefficients,Splines,10, &func_parameters);
    
    // ================ get the avarge of TACs of each segment ==================
    Compute_TACS_Volume(Tvolumes,coefficients, Splines,  tbNum, tbLen, Nvolume);
    
    //free containers with b-spline size to redefine with factors size
    free(coefficients);
    free(Splines);
    
    // Redefine and allocate new continers with the number of segments+1=number of factors
    tbNum = numbSegments;
    Ncoefficients = tbNum * Nvolume;
    
    //declare containers with factors size
    factors = (float *) malloc(size_float * tbLen*tbNum);
    coefficients = (float *) malloc(size_float * Ncoefficients);
    
    func_parameters.basisNum=tbNum;
    func_parameters.Ncoefficients= Ncoefficients;
    func_parameters.Nvolume= Nvolume;
    
    // average tacs using segments
    Compute_Averaged_TACS(Tvolumes,StaticMask,factors,  tbNum,  tbLen, Nvolume);
    // Print curves
    Print_Curves(factors,tbNum,tbLen,"Intermediate TACs Estimated by fitting the b-Splines");
    
    // ============================= Estimate new coefficents for the new curves ======================================================
    
    fprintf(stderr, "=============  Coefficients estimation  =========\n");
    
    // turn regularization functions on now
    func_parameters.SMOOTHNESS_LAMBDA =DEFAULT_SMOOTHNESS_LAMBDA;
    func_parameters.MIX_LAMBDA = DEFAULT_MIX_LAMBDA;
    
    // reconstruct coefficients by MLEM
    Spline_MLEM_recon(coefficients,factors, N_ITERATIONS_MLEM, &func_parameters);
    
    // ============================= refine with FADS ============================================================================
    fprintf(stderr, "=============  FADS refining (MLEM)  =========\n");
    
    //reconstruct coefficients by MLEM
    FADS_MLEM_recon(coefficients,factors,&func_parameters);
    
    // average tacs using segments
    Compute_TACS_Volume(Tvolumes,coefficients, factors,  tbNum, tbLen, Nvolume);
    
    // ============================= Save and Print Results ============================================================================
    // Print final factors
    strcpy(buf, coeffs_fname);
    strcat(buf, "_Final_Factors");
    Print_Curves(factors,tbNum,tbLen,buf);
    
    //Save the final factors into a CSV file
    Save_Curves(factors,  tbNum,  tbLen, buf);
    
    // Save final coefficients
    strcpy(buf, coeffs_fname);
    Save_Volume(coefficients,Ncoefficients,buf,"coff");
    
    // Uncmment this if you want save the final time-dependent volumes.
    //Save_Volume(Tvolumes,tbLen * Nvolume,"Time_dependent_Volume_FADS_3_2rotation","Vol");
    
    // Compute final TACs
    Compute_Averaged_TACS(Tvolumes,StaticMask,factors,  tbNum,  tbLen, Nvolume);
    
    // Print final TACs
    Print_Curves(factors,tbNum,tbLen,"Final TACs");
    
    //Save the final TACs into a CSV file
    strcpy(buf, coeffs_fname);
    strcat(buf, "_Final_TACs");
    Save_Curves(factors,  tbNum,  tbLen, buf);
    
    free(Tvolumes);
    free(StaticMask);
    free(coefficients);
    free(sinogram);
    free(factors);
    return;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------  Static reconstruction  ----------------------------
//---------------------------------------------------------------------------------------------------
// recostructs a static volume with MLEM algrithm + TV smmoothing regularization
void ST_MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, float *volume,int Nx,int Ny, int Nz,float SMOOTHNESS_LAMBDA, char *fname ) {
    int n, niter;
    float *tmp_volume, *tmp_sinogram, *constant_denominator;
    float denominator,no_reg_denominator,reg_denominator,Threshold;
    float *mask, *grad;
    //this are used for tv smoothing
    float NN=0.0;
    
    // ------ create ml_em variables
    mask = (float *) malloc(size_float * Nvolume);
    grad = (float *) malloc(size_float * Nvolume);
    
    tmp_sinogram = (float *) malloc( (size_float * Nsinogram));
    constant_denominator = (float *) malloc( (size_float * Nvolume));
    tmp_volume = (float *) malloc( (size_float * Nvolume));
    
    // --- initial volume assignment: all pixels are one
    for (n = 0; n < Nvolume; n++)
    {
        volume[n] = 1.;
        grad[n]=0.0;
        mask[n]=1.0;
    }
    
    // --- compute  element-by-element inverse of efficiency matrix
    for (n = 0; n < Nsinogram; n++)
        tmp_sinogram[n] = 1.;
    
    RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, constant_denominator, fname);
    
    //  -------- ITERATION LOOP --------
    for (niter = 1; niter <= N_ITERATIONS_MLEM; niter++) {
        fprintf(stderr, "Iteration No [%d] of [%d] ", niter, N_ITERATIONS_MLEM);
        
        // compute the reprojection through the n-1 version of the file into tmp_sinogram
        RayDrFP_SM(Nsinogram, tmp_sinogram, Nvolume, volume, fname);
        
        // divide the sinogram by the tmp_sinogram
        for (n = 0; n < Nsinogram; n++)
        {
            // Divide the sinogram by the estimated sinogram
            if (sinogram[n] == 0.)
                tmp_sinogram[n] = 0.;
            else if (sinogram[n] < 0.)
                fprintf(stderr, "sinogram in MLEM smaller than zero");
            else if (tmp_sinogram[n] > EPSILON)
                tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
            else
                tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
        }
        
        // If smoothness parameter is greater than zero do the smoothness regularization
        if(SMOOTHNESS_LAMBDA>0.0)
        {
            // create a new mask from the current volume by thresholding
            Threshold= otsu( volume,  Nvolume);
            Convert2Mask(mask, volume, Nvolume, Threshold);
            
            // Nearest neighbour constraint
            NN = NN_3D(grad, volume, mask,Nx, Ny, Nz, 100);
            fprintf( stderr,"Smoothness value=:%f: ",NN);
            fprintf( stderr,"TV_LAMBDA=:%f:",SMOOTHNESS_LAMBDA);
        }
        
        // backproject the result into tmp_volume
        RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, tmp_volume, fname);
        
        // multiply by the constant denominator
        for (n = 0; n < Nvolume; n++) {
            no_reg_denominator= (constant_denominator[n]>EPSILON? constant_denominator[n]:EPSILON);
            
            if(SMOOTHNESS_LAMBDA>0.0)
            {
                reg_denominator= (constant_denominator[n] + SMOOTHNESS_LAMBDA * grad[n]);
                
                if (reg_denominator > 0.001)
                {
                    denominator = 1. / reg_denominator;
                }
                else if (no_reg_denominator < 0.001)
                {
                    denominator = 1./no_reg_denominator;
                }
                else
                {
                    denominator = 0.001;
                }
            }
            else
            {
                denominator = 1./no_reg_denominator;
            }
            
            volume[n] *=  denominator * tmp_volume[n];
            
            grad[n] =0.0;
        }
        fprintf(stderr, "\n");
    }
    
    // end: free memory up
    free(mask);
    free(grad);
    free(tmp_sinogram);
    free(tmp_volume);
    free(constant_denominator);
    return;
}

/////////////// Sinogram driven system-matrix based backprojection
void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,
                char *fname) {
    int ns, n, Nchunk;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nvolume) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Nvolume);
        exit(1);
    }
    
    // set volume values to zero
    memset(volume, 0, size_float*Nvolume);
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Main loop
    for (ns = 0; ns < Nsinogram; ns++) {
        // read volume index and compare to expected
        fread(&n, size_int, 1, fid);
        if (n != ns) {
            fprintf(stderr,"Read in sinogram index %d not equal to expected %d\n", n,
                    ns);
            exit(1);
        }
        // fread chunk size and indices and SM chunks
        fread(&Nchunk, size_int, 1, fid);
        if (Nchunk > Nvolume) {
            fprintf(stderr, "Sinogram chunk %d is longer than Nvolume %d\n",Nchunk, Nvolume);
            exit(1);
        }
        fread(Itemp, size_int, Nchunk, fid);
        fread(SMtemp, size_float, Nchunk, fid);
        // is sinogram pixel is non-zero, do loop muptiplication
        if (fabs(sinogram[ns]) > 1.0e-14)
            for (n = 0; n < Nchunk; n++)
            {
                if(Itemp[n] < Nvolume)
                    volume[Itemp[n]] += sinogram[ns] * SMtemp[n];
            }
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

/////////////// Sinogram driven system-matrix based forward projection
void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,
                char *fname) {
    int ns, n, Nchunk;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nvolume) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Nvolume);
        exit(1);
    }
    // set sinogram values to zero
    // memset(sinogram, 0, size_float*Nsinogram); // don't need to do it here, done in the loop
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Main loop
    for (ns = 0; ns < Nsinogram; ns++) {
        // read volume index and compare to expected
        fread(&n, size_int, 1, fid);
        if (n != ns) {
            fprintf(stderr,
                    "Read in sinogram index %d not equal to expected %d\n", n,
                    ns);
            exit(1);
        }
        // fread chunk size and indices and SM chunks
        fread(&Nchunk, size_int, 1, fid);
        if (Nchunk > Nvolume) {
            fprintf(stderr, "ns = %d, nchunk = %d\n", ns, Nchunk);
            fprintf(
                    stderr,
                    "Forward proj: System matrix chunk length %d is longer than Nvolume %d\n",
                    Nchunk, Nvolume);
            exit(1);
        }
        fread(Itemp, size_int, Nchunk, fid);
        fread(SMtemp, size_float, Nchunk, fid);
        sinogram[ns] = 0.;
        // do loop muptiplication
        for (n = 0; n < Nchunk; n++)
            if (fabs(volume[Itemp[n]]) > 1.0e-14)
                if(Itemp[n] < Nvolume)
                    sinogram[ns] += volume[Itemp[n]] * SMtemp[n];
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

//------------------------------------------------------------------------------------------------------------
//---------------------------------------------  MLEM/MAP Dynamic reconstruction  ----------------------------
//------------------------------------------------------------------------------------------------------------
// Spline metod: estimates coefficients of given curves (splines/basis)
void Spline_MLEM_recon(float *coeffVolume,float *basis,int N_iter, struct PARAMETERS *func_parameters)
{
    int n,i, niter;
    int basisLen,basisNum,Nvolume,Nsinogram,No_Heads,N_rotations,Ncoefficients,NumBins,Nx,Ny,Nz;
    float TV=0.0,mix=0.0,denominator,no_reg_denominator,reg_denominator;
    float *sinogram,*tmp_Coeffsvolume, *tmp_sinogram, *constant_denominator,*tmp_constant_denominator, *tmp_basis;
    float *mask,*grad_cSmooth,*grad_cMix ;
    
    sinogram =func_parameters->sinogram;
    basisLen = func_parameters->basisLen;
    basisNum = func_parameters->basisNum;
    Nvolume  = func_parameters->Nvolume;
    Ncoefficients = func_parameters->Ncoefficients;
    NumBins = func_parameters->NumBins;
    Nsinogram = func_parameters->Nsinogram;
    N_rotations = func_parameters->N_rotations;
    No_Heads = func_parameters->No_Heads ;
    Nx = func_parameters->Nx;
    Ny = func_parameters->Ny;
    Nz =func_parameters->Nz;
    
    //-----------------------------------------
    // ------ create ml_em variables:
    mask = (float *) malloc(size_float * Ncoefficients);
    grad_cMix = (float *) malloc(size_float * Ncoefficients);
    grad_cSmooth = (float *) malloc(size_float * Ncoefficients);
    tmp_sinogram = (float *) malloc( (size_float * Nsinogram));
    
    constant_denominator = (float *) malloc((size_float * Ncoefficients));
    tmp_Coeffsvolume = (float *) malloc( (size_float * Ncoefficients));
    tmp_constant_denominator = (float *) malloc( (size_float * Ncoefficients));
    tmp_basis = (float *) malloc( (size_float * basisNum * basisLen));
    memset(tmp_basis,0.0,basisNum * basisLen);
    
    //////////////////////////////////////////////////////////////////////////////////////
    // --- initial volume assignment: all pixels are one
    for (n = 0; n < Ncoefficients; n++)
    {
        coeffVolume[n]=1.0;
        grad_cMix[n]=0.0;
        grad_cSmooth[n]=0.0;
        tmp_Coeffsvolume[n]=1.0;
        mask[n]=0.0;
    }
    
    // --- compute  element-by-element inverse of efficiency matrix
    for (n = 0; n < Nsinogram; n++)
    {
        tmp_sinogram[n]=1.0;
    }
    
    TimeBasisBP_c(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, basisNum, basisLen,
                  tmp_sinogram,constant_denominator,tmp_constant_denominator, basis, func_parameters->sysmat_fname);
    
    //  -------- ITERATION LOOP --------
    for (niter = 1; niter <= N_iter; niter++) {
        
        fprintf(stderr, "Iter No %d of %d ", niter, N_iter);
        
        // Create a normalized coefficients for the mask
        Normailze_coefficients(tmp_Coeffsvolume, tmp_basis,  basisNum,  basisLen,  Nvolume);
        
        // if smoothing is on, start smooting after few iterations
        if(niter>5 && func_parameters->SMOOTHNESS_LAMBDA>0.0)
        {
            TV =0.0;
            for (i = 0; i < basisNum; i++)
                TV += NN_3D(&grad_cSmooth[i * Nvolume],&coeffVolume[i * Nvolume],&mask[i * Nvolume],Nx,Ny,Nz, basisNum+1);
            fprintf( stderr,"Coeff_smooth:%f:",TV);
            fprintf( stderr,"SMOOTH_LAMBDA:%f:",func_parameters->SMOOTHNESS_LAMBDA);
        }
        
        // if mix constraint is on, start it after few iterations
        if(niter>5 && func_parameters->MIX_LAMBDA>0.0)
        {
            CreateMask(mask, func_parameters->StaticMask,tmp_Coeffsvolume, Nvolume,basisNum,func_parameters);
            
            //=====================================================
            // Coefficient mix constraint
            mix = Coeff_Mix(grad_cMix, coeffVolume, Nvolume, basisNum,mask);
            fprintf( stderr,"Coeff_mix:%f:",mix);
            fprintf( stderr,"MIX_LAMBDA:%f:", func_parameters->MIX_LAMBDA);
        }
        fprintf(stderr, "\n");
        
        //===================== Coefficients ======================================================
        // compute the reprojection through the n-1 version of the file into tmp_sinogram
        TimeBasisFP(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, basisNum,
                    basisLen, tmp_sinogram, coeffVolume, basis, func_parameters->sysmat_fname);
        
        
        for (n = 0; n < Nsinogram; n++)
        {
            // divide the sinogram by the tmp_sinogram
            if (sinogram[n] == 0.)
                tmp_sinogram[n] = 0.;
            else if (sinogram[n] < 0.)
                fprintf(stderr, "sinogram in MLEM smaller than zero");
            else if (tmp_sinogram[n] > EPSILON)
                tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
            else
                tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
        }
        
        // backproject the result into tmp_volume
        TimeBasisBP_c(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, basisNum,
                      basisLen, tmp_sinogram, tmp_Coeffsvolume,tmp_constant_denominator, basis, func_parameters->sysmat_fname);
        
        // multiply by the constant denominator
        //===================== Update Coefficients ======================================================
        for (i = 0; i < basisNum; i++)
        {
            for (n = 0; n < Nvolume; n++)
            {
                //Add regularization to the denominator
                reg_denominator= constant_denominator[i*Nvolume+n] + func_parameters->SMOOTHNESS_LAMBDA * grad_cSmooth[i*Nvolume+n]+ func_parameters->MIX_LAMBDA * grad_cMix[i*Nvolume+n];
                
                // denominator without regularization (this to compare with regularized denominator)
                no_reg_denominator= (constant_denominator[i*Nvolume+n]>EPSILON ? constant_denominator[i*Nvolume+n]:EPSILON);
                
                // if denominator after regularization decreasing, make sure it doesn't become negative.
                if (reg_denominator < no_reg_denominator)
                {
                    if (reg_denominator >= 1)
                    {
                        denominator = 1. / reg_denominator;
                    }
                    else
                    {
                        denominator = 1.;
                    }
                }
                // if denominator after regularization increasing, it will always be positve.
                else
                {
                    denominator = 1. / reg_denominator;
                }
                
                
                // multiply by the result
                coeffVolume[i*Nvolume+n] *=  denominator * tmp_Coeffsvolume[i*Nvolume+n];
                
                if(coeffVolume[i*Nvolume+n] < EPSILON)
                    coeffVolume[i*Nvolume+n] = EPSILON;
                
                // zero the regularization matrix for the next iteration
                grad_cMix[i*Nvolume+n]=0.0;
                grad_cSmooth[i*Nvolume+n]=0.0;
                
                // copy the current coeffs to this matrix to be used to create the mask
                tmp_Coeffsvolume[i*Nvolume+n]=coeffVolume[i*Nvolume+n];
            }
        }
        
    }
    
    free(tmp_basis);
    free(tmp_sinogram);
    free(tmp_Coeffsvolume);
    free(constant_denominator);
    free(grad_cMix);
    free(grad_cSmooth);
    free(tmp_constant_denominator);
    free(mask);
        
    return;
}

// FADS method: estimates both coefficients and curves/factors at the same time
void FADS_MLEM_recon(float *coeffVolume,float *factors,struct PARAMETERS *func_parameters)
{
    int i,j,n, niter;
    int FactorsLen,FactorsNum,Nvolume,Nsinogram,No_Heads,N_rotations,Ncoefficients,NumBins,Nx,Ny,Nz;
    float TV=0.0,mix=0.0,TB_smoothness=0.0;
    float denominator,reg_denominator,no_reg_denominator;
    float *sinogram, *tmp_coeffvolume, *tmp_sinogram,*tmp_basis, *constant_denominator_c, *constant_denominator_f;
    float *Dynamic_mask,*grad_cSmooth,*grad_cMix ,*grad_TB_Smooth;
    
    sinogram =func_parameters->sinogram;
    FactorsLen = func_parameters->basisLen;
    FactorsNum = func_parameters->basisNum;
    Nvolume  = func_parameters->Nvolume;
    Ncoefficients = func_parameters->Ncoefficients;
    NumBins = func_parameters->NumBins;
    Nsinogram = func_parameters->Nsinogram;
    N_rotations= func_parameters->N_rotations;
    No_Heads = func_parameters->No_Heads ;
    Nx = func_parameters->Nx;
    Ny = func_parameters->Ny;
    Nz =func_parameters->Nz;
    
    //-----------------------------------------
    // ------ create ml_em variables:
    tmp_sinogram = (float *) malloc( (size_float * Nsinogram));
    constant_denominator_c = (float *) malloc((size_float * Ncoefficients));
    tmp_coeffvolume = (float *) malloc((size_float * Ncoefficients));
    Dynamic_mask = (float *) malloc(size_float * Ncoefficients);
    grad_cMix = (float *) malloc(size_float * Ncoefficients);
    grad_cSmooth = (float *) malloc(size_float * Ncoefficients);
    grad_TB_Smooth = (float *) malloc(size_float * FactorsNum*FactorsLen);
    tmp_basis = (float *) malloc(size_float * FactorsLen* FactorsNum);
    constant_denominator_f = (float *) malloc(size_float * FactorsLen* FactorsNum);
    //////////////////////////////////////////////////////////////////////////////////////
    for (n = 0; n < Ncoefficients; n++)
    {
        grad_cMix[n]=0.0;
        grad_cSmooth[n]=0.0;
    }
    memset(grad_TB_Smooth, 0, size_float*FactorsNum*FactorsLen);
    
    //  -------- ITERATION LOOP --------
    for (niter = 1; niter <= N_ITERATIONS_MLEM; niter++) {
                fprintf(stderr, "Iter No %d of %d :", niter, N_ITERATIONS_MLEM);
        
        // Normalize the coefficients
        Normailze_coefficients(coeffVolume, factors,  FactorsNum,  FactorsLen,  Nvolume);
        
        //=====================================================
        // Create a maske from current estimated coefficients
        
        CreateMask(Dynamic_mask,func_parameters->StaticMask, coeffVolume, Nvolume,FactorsNum,func_parameters);
        
        //=======================================
        // Nearest neighbour constraint
        TV=0;
        for (i = 0; i < FactorsNum; i++) {
            TV += NN_3D(&grad_cSmooth[i * Nvolume], &coeffVolume[i * Nvolume],&Dynamic_mask[i * Nvolume],Nx,Ny,Nz, FactorsNum+1);
        }
                        fprintf( stderr,"Coeff_smooth:%f:",TV);
                        fprintf( stderr,"SMOOTH_LAMBDA:%f:",func_parameters->SMOOTHNESS_LAMBDA);
        
        //=====================================================
        // Coefficient mix constraint
        mix = Coeff_Mix(grad_cMix, coeffVolume, Nvolume, FactorsNum,Dynamic_mask);
        
                fprintf( stderr,"Coeff_Mix:%f:",mix);
                fprintf( stderr,"MIX_LAMBDA:%f:", func_parameters->MIX_LAMBDA);
        
        //=======================================
        //factors smoothness constraint
        TB_smoothness = TB_Smooth(grad_TB_Smooth,factors,  FactorsLen,  FactorsNum);
                        fprintf( stderr,"TB_smooth:%f:",TB_smoothness);
                        fprintf( stderr,"TB_SMOOTH_LAMBDA:%f:\n",func_parameters->TB_SMOOTHNESS_LAMBDA);
        
        
        // compute the reprojection through the n-1 version of the file into tmp_sinogram
        TimeBasisFP(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, FactorsNum,
                    FactorsLen, tmp_sinogram, coeffVolume, factors, func_parameters->sysmat_fname);
        
        for (n = 0; n < Nsinogram; n++)
        {
            // divide the sinogram by the tmp_sinogram
            if (sinogram[n] == 0.)
                tmp_sinogram[n] = 0.;
            else if (sinogram[n] < 0.)
                fprintf(stderr, "sinogram in MLEM smaller than zero");
            else if (tmp_sinogram[n] > EPSILON)
                tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
            else
                tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
        }
        
        // backproject the result into tmp_volume
        TimeBasisBP_c(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, FactorsNum,
                      FactorsLen, tmp_sinogram, tmp_coeffvolume, constant_denominator_c, factors, func_parameters->sysmat_fname);
        
        //===================== Update Coefficients ======================================================
        for (i = 0; i < Ncoefficients; i++)
        {
            //Add regularization to the denominator
            reg_denominator= constant_denominator_c[i]+func_parameters->SMOOTHNESS_LAMBDA*grad_cSmooth[i] +func_parameters->MIX_LAMBDA * grad_cMix[i];
            no_reg_denominator= (constant_denominator_c[i]>EPSILON ? constant_denominator_c[i]:EPSILON);
            
            //  apply the regularization
            if (reg_denominator < no_reg_denominator)
            {
                if (reg_denominator >= 1)
                {
                    denominator = 1. / reg_denominator;
                }
                else
                {
                    denominator = 1.;
                }
            }
            else
            {
                denominator = 1. / reg_denominator;
            }
            
            
            // update coefficients
            coeffVolume[i] *=  denominator * tmp_coeffvolume[i];
            if(coeffVolume[i] < EPSILON)
                coeffVolume[i] = EPSILON;
            
            // zero the regularization matrix for the next iteration
            grad_cMix[i]=0.0;
            grad_cSmooth[i]=0.0;
        }
        
        // compute the reprojection through the n-1 version of the file into tmp_sinogram
        TimeBasisFP(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, FactorsNum,
                    FactorsLen, tmp_sinogram, coeffVolume, factors, func_parameters->sysmat_fname);
        
        // divide the sinogram by the tmp_sinogram
        for (n = 0; n < Nsinogram; n++)
        {
            // divide the sinogram by the tmp_sinogram
            if (sinogram[n] == 0.)
                tmp_sinogram[n] = 0.;
            else if (sinogram[n] < 0.)
                fprintf(stderr, "sinogram in MLEM smaller than zero");
            else if (tmp_sinogram[n] > EPSILON)
                tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
            else
                tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
        }
        
        
        // backproject the result into tmp_volume
        TimeBasisBP_f(Nsinogram,No_Heads, N_rotations, Ncoefficients, NumBins, FactorsNum,
                      FactorsLen, tmp_sinogram, coeffVolume, constant_denominator_f, tmp_basis, func_parameters->sysmat_fname);
        
        //===================== Update factors ======================================================
        for (i = 0; i < (FactorsNum); i++)
        {
            for (j = 0; j < FactorsLen; j++)
            {
                //Add regularization to the denominator
                reg_denominator= (constant_denominator_f[j+i*FactorsLen] + (func_parameters->TB_SMOOTHNESS_LAMBDA * grad_TB_Smooth[j+i*FactorsLen]));
                no_reg_denominator= (constant_denominator_f[j+i*FactorsLen]>EPSILON ? constant_denominator_f[j+i*FactorsLen]:EPSILON);
                
                //  apply the regularization
                if (reg_denominator < no_reg_denominator)
                {
                    if (reg_denominator >= 1)
                    {
                        denominator = 1. / reg_denominator;
                    }
                    else
                    {
                        denominator = 1.;
                    }
                }
                else
                {
                    denominator = 1. / reg_denominator;
                }
                
                
                factors[j+i*FactorsLen] *=  denominator * tmp_basis[j+i*FactorsLen];
                
                if(factors[j+i*FactorsLen] <= EPSILON) factors[j+i*FactorsLen] = EPSILON;
                
                // zero the regularization matrix for the next iteration
                grad_TB_Smooth[j+i*FactorsLen]=0.0;
            }
            
            // Smooth curve by fft
            Smooth_curve( &factors[i*FactorsLen],  FactorsLen);
            
            
        }

    }
    
    free(grad_TB_Smooth);
    free(grad_cMix);
    free(grad_cSmooth);
    free(Dynamic_mask);
    free(constant_denominator_f);
    free(tmp_basis);
    free(tmp_sinogram);
    free(tmp_coeffvolume);
    free(constant_denominator_c);
    
    return;
}

/* back projection using a system-matrix, modified by a set of time basis functions */
void TimeBasisBP_c(int Nsinogram,int No_Heads, int N_rotations, int Ncoefficients, int NumBins,
                   int basisNum, int basisLen, float *sinogram, float *coeffVolume, float *denominator_c,
                   float *basis, char *sysmat_fname) {
    int  b, j,r,np,h,n, Nchunk, Nvolume,Nsino_1h_1r, Number_Projections_1Rotation_1head;
    int coeffs_index, sinogram_index, basis_index;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n * N_rotations != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // size of sinogram of one head and one rotation
    Nsino_1h_1r = n/No_Heads;
    
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if ((n * basisNum) != Ncoefficients) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Ncoefficients);
        exit(1);
    }
    
    Nvolume = n;
    // set volume values to zero
    for (j = 0; j < Ncoefficients; j++) {
        coeffVolume[j]=0.0;
        denominator_c[j]=0.0;
    }
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Number of projections in one Rotation and one head
    Number_Projections_1Rotation_1head =basisLen/N_rotations;
    
    // Main loop
    for (h = 0; h < No_Heads; h++) {
        for (np = 0; np < Number_Projections_1Rotation_1head; np++) {
            for (b = 0; b < NumBins; b++) {
                
                //calcalute sinogram index
                sinogram_index = b+np*NumBins+ h*NumBins*Number_Projections_1Rotation_1head;
                
                // read sinogram index and compare to expected
                fread(&n, size_int, 1, fid);
                if (n != sinogram_index) {
                    fprintf(stderr,
                            "Read in sinogram index %d not equal to expected %d\n", n,sinogram_index);
                    exit(1);
                }
                // fread chunk size and indices and SM chunks
                fread(&Nchunk, size_int, 1, fid);
                if (Nchunk > Nvolume) {
                    fprintf(stderr, "Sinogram chunk %d is longer than Nvolume %d\n",
                            Nchunk, Nvolume);
                    exit(1);
                }
                fread(Itemp, size_int, Nchunk, fid);
                fread(SMtemp, size_float, Nchunk, fid);
                
                for (r = 0; r < N_rotations; r++) {
                    for (n = 0; n < Nchunk; n++) {
                        for (j = 0; j < basisNum; j++) {
                            
                            //calcalute sinogram index
                            sinogram_index = b+np*NumBins + r*Nsino_1h_1r + h*NumBins*basisLen;
                            
                            //calcalute coefficients index
                            coeffs_index= Itemp[n] + j * Nvolume;
                            
                            //calcalute basis index
                            basis_index= np + r*Number_Projections_1Rotation_1head + j*basisLen ;
                            
                            if( basis[basis_index]>EPSILON)
                            {
                                if(fabs(sinogram[sinogram_index])>EPSILON)
                                    coeffVolume[coeffs_index] += sinogram[sinogram_index] * SMtemp[n] * basis[basis_index];
                                
                                denominator_c[coeffs_index] +=  SMtemp[n] * basis[basis_index];
                            }
                        }
                    }
                }
            }
        }
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

/* back projection using a system-matrix, modified by a set of time basis functions */
void TimeBasisBP_f(int Nsinogram,int No_Heads,int N_rotations, int Ncoefficients, int NumBins,
                   int basisNum, int basisLen, float *sinogram, float *coeffVolume,
                   float *denominator_f, float *basis, char *sysmat_fname) {
    int i,j, k,l, m,n,h, Nchunk, Nvolume,Nsino_1h_1r, Number_Projections_1Rotation_1head;
    int coeffs_index, sinogram_index, basis_index;
    float sum;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n * N_rotations != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // size of sinogram of one head and one rotation
    Nsino_1h_1r = n/No_Heads;
    
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if ((n * basisNum) != Ncoefficients) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Ncoefficients);
        exit(1);
    }
    
    Nvolume = n;
    // set values to zero
    
    for (k = 0; k < basisNum*basisLen; k++) {
        basis[k]=0.0;
        denominator_f[k]=0.0;
    }
    
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Number of projections in one Rotation and one head
    Number_Projections_1Rotation_1head =basisLen/N_rotations;
    
    // Main loop
    for (h = 0; h < No_Heads; h++) {
        for (l = 0; l < Number_Projections_1Rotation_1head; l++) {
            for (i = 0; i < NumBins; i++) {
                
                //calcalute sinogram index
                sinogram_index = i+l*NumBins+ h*NumBins*Number_Projections_1Rotation_1head;
                
                // read sinogram index and compare to expected
                fread(&n, size_int, 1, fid);
                if (n != sinogram_index) {
                    fprintf(stderr,
                            "Read in sinogram index %d not equal to expected %d\n", n,sinogram_index);
                    exit(1);
                }
                // fread chunk size and indices and SM chunks
                fread(&Nchunk, size_int, 1, fid);
                if (Nchunk > Nvolume) {
                    fprintf(stderr, "Sinogram chunk %d is longer than Nvolume %d\n",
                            Nchunk, Nvolume);
                    exit(1);
                }
                fread(Itemp, size_int, Nchunk, fid);
                fread(SMtemp, size_float, Nchunk, fid);
                
                // is sinogram pixel is non-zero, do loop muptiplication
                for (k = 0; k < (basisNum); k++)
                {
                    sum = 0.0;
                    for (j = 0; j < Nchunk; j++) {
                        
                        //calcalute coefficients index
                        coeffs_index= Itemp[j] + k * Nvolume;
                        
                        if(coeffVolume[coeffs_index]>EPSILON)
                            sum += (SMtemp[j] * coeffVolume[coeffs_index]);
                    }
                    
                    for (m = 0; m < N_rotations; m++) {
                        
                        //calcalute sinogram index
                        sinogram_index = i+l*NumBins + m*Nsino_1h_1r + h*NumBins*basisLen;
                        
                        //calcalute basis index
                        basis_index= k*basisLen + l + m*Number_Projections_1Rotation_1head;
                        
                        if((sinogram[sinogram_index])>EPSILON)
                            basis[basis_index] += sinogram[sinogram_index] * sum;
                        
                        denominator_f[basis_index] += sum;
                    }
                }
            }
        }
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

/* Forward projection using a system-matrix, modified by a set of time basis functions */
void TimeBasisFP(int Nsinogram,int No_Heads, int N_rotations, int Ncoefficients, int NumBins,
                 int basisNum, int basisLen, float *sinogram, float *coeffVolume,
                 float *basis, char *sysmat_fname) {
    int i, j, k,l, n, h, Nchunk, Nvolume,Nsino_1h_1r, Number_Projections_1Rotation_1head;
    int coeffs_index, sinogram_index, basis_index;
    float *SMtemp;
    int *Itemp;
    float sum;
    FILE *fid;
    
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", sysmat_fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n * N_rotations != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // size of sinogram of one head and one rotation
    Nsino_1h_1r = n/No_Heads;
    
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if ((n * basisNum) != Ncoefficients) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Ncoefficients);
        exit(1);
    }
    
    Nvolume = n;
    
    // set sinogram values to zero
    memset(sinogram, 0, size_float*Nsinogram);
    
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Number of projections in one Rotation and one head
    Number_Projections_1Rotation_1head =basisLen/N_rotations;
    
    // Main loop
    for (h = 0; h < No_Heads; h++) {
        for (l = 0; l < Number_Projections_1Rotation_1head; l++) {
            for (i = 0; i < NumBins; i++) {
                
                //calcalute sinogram index
                sinogram_index = i + l * NumBins + h*NumBins*Number_Projections_1Rotation_1head;
                
                // read sinogram index and compare to expected
                fread(&n, size_int, 1, fid);
                if (n != sinogram_index) {
                    fprintf(stderr,
                            "Read in sinogram index %d not equal to expected %d\n",
                            n, sinogram_index);
                    exit(1);
                }
                // fread chunk size and indices and SM chunks
                fread(&Nchunk, size_int, 1, fid);
                if (Nchunk > Nvolume) {
                    fprintf(stderr, "ns = %d, nchunk = %d\n", l, Nchunk);
                    fprintf(
                            stderr,
                            "Forward proj: System matrix chunk length %d is longer than Nvolume %d\n",
                            Nchunk, Nvolume);
                    exit(1);
                }
                
                fread(Itemp, size_int, Nchunk, fid);
                fread(SMtemp, size_float, Nchunk, fid);
                
                // do loop muptiplication
                for (j = 0; j < basisNum; j++)
                {
                    sum = 0.0;
                    for (n = 0; n < Nchunk; n++)
                    {
                        
                        //calcalute coefficients index
                        coeffs_index= Itemp[n] + j * Nvolume;
                        
                        if(coeffVolume[coeffs_index]>EPSILON)
                            sum += coeffVolume[coeffs_index] * SMtemp[n];
                    }
                    
                    for (k = 0; k < N_rotations; k++)
                    {
                        //calcalute sinogram index
                        sinogram_index = i + l * NumBins  + k*Nsino_1h_1r + h*NumBins*basisLen;
                        
                        //calcalute basis index
                        basis_index= j*basisLen+ l + k*Number_Projections_1Rotation_1head;
                        
                        if (basis[basis_index]>EPSILON)
                            sinogram[sinogram_index] += sum * basis[basis_index];
                    }
                }
            }
        }
    }
    
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

//---------------------------------------------------------------------------------------------------
//--------------------------------------------- Regularization functions ----------------------------
//---------------------------------------------------------------------------------------------------

// Nearest neighbours smoothing constraint
float NN_3D(float *NN_grad, float *volume,float *mask, int Nx, int Ny, int Nz, int omega)
{
    int x, y, z,i,j,k,xx,yy,zz, vox, N_vox;
    float diff,abs_diff;
    float NNV=0.0 ;
    
    for(z=0; z<Nz; z++)
        for(y=0; y<Ny; y++)
            for(x=0; x<Nx; x++)
            {
                // find the current voxel position in the 1d array
                vox = voxel_pos(Nx,Ny, x,y,z);
                
                // initial value
                NN_grad[vox]=0.0;
                
                // 27 nearest voxels.
                for(i=0; i<3; i++) // i is index the silce (before , same, after the voxel's slice)
                {
                    for(j=0; j<3; j++) // j is index of the current clice in x direction.
                    {
                        for(k=0; k<3; k++) // k is index of the current clice in y direction.
                        {
                            // get the current nearest voxel coordinates
                            xx= x+(j-1);
                            yy= y+(k-1);
                            zz= z+(i-1);
                            
                            // skip if the nearest voxel is out of range (the current voxel could be at one outer slices of the volume)
                            if( xx<0 || xx>(Nx-1) || yy<0 || yy>(Ny-1) || zz<0 || zz>(Nz-1))
                                continue;
                            
                            // find the current nearest voxel position in the 1d array
                            N_vox = voxel_pos(Nx,Ny, xx,yy,zz);
                            
                            // if this voxel is not voxel that has mix ( if it is a mix voxel then it will penalized for mix not for smoothness)
                            // and this voxel located in the same tissue as its neighbour, then apply smoothing penalty
                            if(mask[vox] != omega && mask[vox] == mask[N_vox])
                            {
                                diff = (volume[vox] - volume[N_vox]);
                                abs_diff = fabs(diff);
                                
                                // if there is a vaiation calculate the function value and its derivative
                                if(abs_diff>EPSILON)
                                {
                                    NNV += fabs(diff);
                                    NN_grad[vox] +=  (diff/abs_diff);
                                }
                            }
                        }
                    }
                }
            }
    
    // check if overflow occured
    if(isnan(NNV))
    {
        fprintf(stderr,"\n \n 2- Error: overflow occured in coefficients smoothness constraint computation !!! \n");
        exit(1);
    }
    return NNV;
}

float Coeff_Mix(float *MixGrad, float *xt,int Nvolume, int basisNum,float *mask)
{
    int i,k,l;
    float MixFunction=0.0;
    float Dot_product_prime;
    
    
    // coefficients mix constraint
    Dot_product_prime = 0.;
    for (l = 0; l < Nvolume; l++)
    {
        // calculate the dot product
        for (k = 0; k < basisNum; k++)
        {
            MixGrad[l + k * Nvolume]=0.0;
            
            for (i = k+1; i < basisNum; i++)
            {
                MixFunction +=  (xt[l + i * Nvolume] *  xt[l + k * Nvolume]);
            }
            
            // this to calculate the derivayive
            Dot_product_prime += xt[l + k * Nvolume];
            
        }
        
        // compute the penlty and derivative of the function with respect to each coefficent that has mix
        
        for (k = 0; k < basisNum; k++)
        {
            for (i = 0; i < basisNum; i++)
            {
                if(k==i) continue;
                
                if (mask[l + k * Nvolume]==basisNum+1)
                    MixGrad[l + k * Nvolume] +=   (xt[l + i * Nvolume]);
            }
        }
    }
    
    //check if there is overflow in calculation
    if(isnan(MixFunction))
    {
        fprintf(stderr,"\n \n Error: overflow occured in coefficients mix constraint computation !!! \n");
        exit(1);
    }
    
    return MixFunction;
}

// Factor Smoothness constraint
float TB_Smooth(float *TBSmooth_grad,float *basis, int basisLen, int basisNum)
{
    int i, j;
    float func, func_prime;
    float Total_TBSmooth_Func=0.0,f0,f1,f2;
    
    for (i = 0; i < basisNum; i++)
    {
        for (j = 0; j < basisLen; j++)
        {
            //first put zero
            TBSmooth_grad[j+i*basisLen] =0.0;
            
            // previous point on the curve (if we are at the first point there is no previous point , so it equal to the current)
            f0 = (j==0)? basis[j+i*basisLen]:basis[j-1+i*basisLen];
            
            // Current point on the curve
            f1 = basis[j+i*basisLen];
            
            // Next point on the curve (if we are at the end point there is no next point, so it equal to the current)
            f2 =  (j==basisLen-1)? basis[j+i*basisLen]:basis[j+1+i*basisLen];
            
            // calculate the total penalty
            func = f2*f1 +f0*f1 - f1*f1;
            func_prime = f2+f0-2.0*f1;
            // check if the total penalty is significant
            if (fabs(func)>0.0)
            {
                // add it to the smoothness function
                Total_TBSmooth_Func += fabs(func);
                // Calculate the derivative of the smoothness function
                TBSmooth_grad[j+i*basisLen] += (func_prime * func)/fabs(func);
            }
        }
    }
    
    return Total_TBSmooth_Func;
}

//---------------------------------------------------------------------------------------------------
//--------------------------------------------- mask creation functions -----------------------------
//---------------------------------------------------------------------------------------------------


// Static Masking: creates a mask out of coefficients and segmented static volume
void CreateMask(float *mask, float *Static_Mask, float *coefficients,int Nvolume,int basisNum,struct PARAMETERS *func_parameters)
{
    int k,l;
    float max;
    int maxIndex;
    float sum,Threshold;
    
    float *convolved_coeffs;
    convolved_coeffs = (float *) malloc(size_float * Nvolume*basisNum);
    
    //Create initial mask from the coefficients volume by tresholding
    for (k = 0; k < basisNum; k++) {
        Threshold= otsu( &coefficients[k * Nvolume],  Nvolume);
        Convert2Mask(&mask[k * Nvolume], &coefficients[k * Nvolume], Nvolume, Threshold);
    }
    
    
    //Average each coefficient with its neighbours
    for (k = 0; k < basisNum; k++) {
        convolve_coeffs(&convolved_coeffs[k * Nvolume], &coefficients[k * Nvolume], &Static_Mask[k * Nvolume],Nvolume,func_parameters);
    }
    
    for (l = 0; l < Nvolume; l++)
    {
        //these two variables are to find the index of the coefficient that has the maximum neighboring average
        max=0.;
        maxIndex= basisNum+3; // initial value is an index out of range
        
        //check if there is mix
        sum=0.;
        for (k = 0; k < basisNum; k++)
        {
            
            if (max<convolved_coeffs[l + k * Nvolume])
            {
                max=convolved_coeffs[l + k * Nvolume];
                maxIndex= k;
            }
            
            //this to check if the mask has more than one coeffcient with value 1 (meaning there is a mix)
            sum += mask[l + k * Nvolume];
        }
        
        // if there is a mix
        if(sum>1)
            for (k = 0; k < basisNum; k++)
            {
                if (maxIndex == k)
                    mask[l + k * Nvolume]=1;
                else
                    mask[l + k * Nvolume]=basisNum+1;  // marks this coefficient to be penalized
            }
    }
    
    free(convolved_coeffs);
    return;
}

//converts a volume into a binary mask according to a given threshhold
void Convert2Mask(float *mask, float *volume,int Nvolume,float Threshold)
{
    int i;
    for (i = 0; i < Nvolume; i++) {
        if(volume[i] >Threshold)
            mask[i]=1.0;
        else
            mask[i]=0.0;
    }
    return;
}

// Averages each coefficient with its neighbours
void convolve_coeffs(float *Ave_NN, float *coefficients,float *mask, int Nvolume,struct PARAMETERS *func_parameters)
{
    int i,j,l,m;
    int x,y,z,xx,yy,zz,vox,N_vox;
    int Nx = func_parameters->Nx;
    int Ny = func_parameters->Ny;
    int Nz = func_parameters->Nz;
    float Nsum;
    int count;
    
    
    for (l = 0; l < Nvolume; l++)
    {
        //index of the coefficient
        vox = l;
        
        // coordinates of the coefficient
        z= l/(Nx*Ny);
        y=(l- z * (Nx*Ny))/Nx;
        x=(l- z * (Nx*Ny))-y*Nx;
        
        // sum of 27 nearest coefficient.
        Nsum=0.0;
        count=0;
        for(i=0; i<3; i++) // i is index the silce (before , same, after the voxel's slice)
        {
            for(j=0; j<3; j++) // j is index of the current clice in x direction.
            {
                for(m=0; m<3; m++) // k is index of the current clice in y direction.
                {
                    
                    // get the current nearest voxel (coefficient) coordinates
                    xx= x+(j-1);
                    yy= y+(m-1);
                    zz= z+(i-1);
                    
                    // skip if the nearest voxel is out of range (the current voxel could be at the one of the outer slices of the volume)
                    if( xx<0 || xx>(Nx-1) || yy<0 || yy>(Ny-1) || zz<0 || zz>(Nz-1))
                        continue;
                    
                    // find the current nearest voxel position in the 1d array
                    N_vox = voxel_pos(Nx,Ny, xx,yy,zz);
                    
                    if(mask[vox] == mask[N_vox])
                    {
                        // if bot the currnet voxel and the nearst voxel in the same region, calculate the deference
                        Nsum += coefficients[N_vox];
                        count++;
                    }
                }
            }
        }
        if(count>0)
            Ave_NN[l]= Nsum/(count);
        else
            Ave_NN[l]=0.;
    }
    return;
}


// finds a threshold of an image using otsu method
float otsu(float *volume, int size)
{
    // NOTE: Creation of histogram[256] not shown
    int   histogram[256];
    
    float  w = 0;                // first order cumulative
    float  u = 0;                // second order cumulative
    float  uT = 0;               // total mean level
    
    int    k = 255;              // maximum histogram index
    float    threshold = 0;        // optimal threshold value
    
    float  histNormalized[256];  // normalized histogram values
    
    float  work1, work2;		// working variables
    double work3 = 0.0;
    float cmax=0;
    int i,j,G;
    
    for(k=0; k<256; k++)  histogram[k]=0;
    
    for ( i = 0; i < size; ++i)
    {
        if (cmax < volume[i])
            cmax = volume[i];
    }
    
    for(i=0; i<size; i++)
    {
        G=(volume[i]/cmax) * 256;
        histogram[G]++;
    }
    
    // Create normalised histogram values
    // (size=image width * image height)
    for ( j=1; j<=k; j++)
        histNormalized[j-1] = histogram[j-1]/(float)size;
    
    
    // Calculate total mean level
    for ( j=1; j<=k; j++)
        uT+=(j*histNormalized[j-1]);
    
    
    // Find optimal threshold value
    for ( j=1; j<k; j++) {
        w+=histNormalized[j-1];
        u+=(j*histNormalized[j-1]);
        work1 = (uT * w - u);
        work2 = (work1 * work1) / ( w * (1.0f-w) );
        if (work2>work3) work3=work2;
    }
    
    // Convert the final value to an integer
    threshold = sqrt(work3)/256.0;
    
    return threshold;
}

//---------------------------------------------------------------------------------------------------
//--------------------------------------------- smooth curve by fft functions -----------------------
//---------------------------------------------------------------------------------------------------

void Smooth_curve( float *curve, int basisLen)
{
    int i,j, count,ave_window=4;
    float sum_a,sum_b;
    int length;
    float *function;
    
    length = basisLen + basisLen/4;
    function = (float *) malloc(size_float *length);
    
    
    for ( j = 0; j < basisLen; ++j)
    {
        function[j]= curve[j];
    }
    
    for ( j = basisLen; j < length; ++j)
    {
        function[j]= curve[j- basisLen/4];
    }
    
    float *a,*b; // real and imaginary
    
    a = (float *) malloc(size_float * length/2+1);
    b = (float *) malloc(size_float * length/2+1);
    
    // compute fft of the curve
    forwardDFT_1D(function, length, a, b);  // forward DFT
    
    // smooth fft curve with a moving window
    for ( j = 5; j < (length/2+1); ++j)
    {
        sum_a= a[j];
        sum_b= b[j];
        count=1;
        for ( i = 0; i < ave_window/2; ++i)
        {
            if((j+i)<(basisLen/2+1))
            {
                sum_a += a[j+i];
                sum_b += b[j+i];
                count++;
                
            }
            
            if((j-i)>1)
            {
                sum_a += a[j-i];
                sum_b += b[j-i];
                count++;
            }
        }
        
        a[j]= sum_a/count;
        b[j]= sum_b/count;
    }
    
    // invert back the smoothed fft curve
    inverseDFT_1D(a, b, length, function); // inverse DFT
    
    for ( j = 0; j < basisLen; ++j)
    {
        if(function[j] <= EPSILON)
            curve[j] = EPSILON;
        else
            curve[j] = function[j];
    }
    
    free(function);
    free(a);
    free(b);
    return;
}

void forwardDFT_1D(const float *data, const int N, float *real, float *imaginary)
{
    int i,j;
    for ( i = 0; i <= N / 2; ++i) {
        real[i] = imaginary[i] = 0;
        for ( j = 0; j < N; ++j) {
            real[i] += data[j] * cos(2 * M_PI / N * i * j);
            imaginary[i] += data[j] * sin(2 * M_PI / N * i * j);
        }
        // normalization
        real[i] *= (i == 0 || i == N / 2) ? 1. / N : 2. / N;
        imaginary[i] *= 2. / N;
    }
}

void inverseDFT_1D(const float *real, const float *imaginary, const int N, float *data)
{
    int i,j;
    for ( j = 0; j < N; ++j) {
        data[j] = real[0];
        for ( i = 1; i <= N / 2; ++i) {
            data[j] += real[i] * cos(2 * M_PI / N * i * j) + imaginary[i] * sin(2 * M_PI / N * i * j);
        }
    }
}

//---------------------------------------------------------------------------------------------------
//--------------------------------------------- help functions --------------------------------------
//---------------------------------------------------------------------------------------------------

// reads volume or sinogram from file to container volume in memory
void Read_volume(float *Data, int read_size,int Data_Type,  char *fname)
{
    FILE *fid;
    int n;
    unsigned short tmp;
    
    // open the file for reading
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr,
                "Could not open volume data file \"%s\" for reading\n",
                fname);
        exit(1);
    }
    
    // read the file
    if(Data_Type == 1) // 8-bit data
    {
        for (n = 0; n < read_size; n++) {
            fread(&tmp, 2, 1, fid); // 2 = sizeof(unsigned short) == 16 bit
            Data[n] = tmp;
        }
    }
    else if (Data_Type == 2) // 32-bit
    {
        fread(Data, size_float, read_size, fid);
    }
    else
    {
        fprintf(stderr,"Can't read data. Data type must either 8-bit (1) or 32-bit (2)\n");
        exit(1);
    }
    fclose(fid);
    
    return;
}

// Saves volume or singram to HDD
void Save_Volume(float *volume, int Save_size,char Save_Name[128],char Save_Extension[128])
{
    FILE *fin;
    char buf[128];
    
    // construct the save name
    strcpy(buf, Save_Name);
    strcat(buf, ".");
    strcat(buf, Save_Extension);
    
    //open file for writing
    if ((fin = fopen(buf, "wb")) == NULL) {
        fprintf(stderr, "Could not open data file %s for writing\n",buf);
        exit(1);
    }
    // write the data
    fwrite(volume, size_float, Save_size, fin);
    
    //close the file
    fclose(fin);
    
    return;
}

// Reads temporal curves (B-splines, Factors, or any curves)
// Format must be:
// first line of the file: number of curves
// second line of the file: length of the curves
// the rest of file contains columns separated by a space or a tab. Each column represents one curve.
void Read_Curves(float *basis, int tbNum,int tbLen,  char *basis_fname)
{
    
    FILE *fid;
    int i,j;
    float temp;
    
    // read the contents of the time basis file into an array
    if ((fid = fopen(basis_fname, "r")) == NULL) {
        fprintf(stderr, "Could not open basis file %s for reading.\n",
                basis_fname);
        exit(1);
    }
    // first read the number and length of time basis
    fscanf(fid, "%d\n", &tbNum);
    fscanf(fid, "%d\n", &tbLen);
    
    
    for(i=0;i<tbLen;i++)
    {
        for(j=0;j<tbNum;j++)
        {
            fscanf(fid, "%f ", &temp);
            basis[j*tbLen+i]=temp;
        }
    }
    fclose(fid);
}

// Saves the curves into a CSV file
void Save_Curves(float *basis, int tbNum, int tbLen, char Save_Name[128] )
{
    FILE *fp;
    char buf[128];
    int i,j;
    
    // construct the save name
    strcpy(buf, Save_Name);
    strcat(buf, ".csv");
    
    //open file for writing
    if ((fp = fopen(buf, "w+")) == NULL) {
        fprintf(stderr, "Could not open text file %s for writing\n",buf);
        exit(1);
    }
    
    // print the curve in a table format
    for (i = 0; i < tbLen; i++) {
        for (j = 0; j < tbNum; j++) {
            
            fprintf(fp, "%f,",basis[j*tbLen+i]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    return;
}

// Prints the curves on the screen
void Print_Curves(float *basis, int tbNum, int tbLen, char title[128] )
{
    int i,j;
    
    // print the header
    fprintf(stderr, "\n============= %s  =========\n",title);
    
    // print the curve in a table format
    for (i = 0; i < tbLen; i++) {
        for (j = 0; j < tbNum; j++) {
            
            fprintf(stderr, "%14.14f	",basis[j*tbLen+i]);
        }
        fprintf(stderr, "\n");
    }
    
    fprintf(stderr, "\n");
    return;
}
// computes the time-dependent volumes (TACS) by multiplying the coefficients by the curves (Cxf)
void Compute_TACS_Volume(float *Tvolumes,float *volume, float *basis, int tbNum, int tbLen,int Nvolume)
{
    int i,j,k;
    float sum;
    
    for (j = 0; j < tbLen; j++) {
        for (i = 0; i < Nvolume; i++) {
            sum=0.0;
            for (k = 0; k< tbNum; k++) {
                if(volume[i + k * Nvolume]>EPSILON && basis[k*tbLen+j]>EPSILON)
                    sum += volume[i + k * Nvolume] *  basis[k*tbLen+j];
            }
            Tvolumes[i + j * Nvolume] = sum ;
        }
    }
}

// averages the given TACS according to given segements
void Compute_Averaged_TACS(float *Tvolumes,float *segmentedVolume, float *basis, int tbNum, int tbLen,int Nvolume)
{
    int i,j,k;
    int Segment_count;
    float sum;
    for (j = 0; j < tbLen; j++) {
        for (k = 0; k< tbNum; k++) {
            Segment_count=0;
            basis[k*tbLen+j]= sum = 0.;
            //sum up this segment and cout how many voxels
            for (i = 0; i < Nvolume; i++) {
                if (segmentedVolume[i+k*Nvolume]==1.0)
                {
                    sum +=  Tvolumes[i + j * Nvolume];
                    Segment_count++;
                }
            }
            
            // average the summed values
            if(Segment_count>0)
            {
                basis[k*tbLen+j] = sum /Segment_count;
                if (basis[k*tbLen+j]<= EPSILON) basis[k*tbLen+j]=EPSILON;
            }
            else
            {
                basis[k*tbLen+j] =EPSILON;
            }
        }
    }
    
    return;
}

//Normalizes the coefficients
void Normailze_coefficients(float *coefficients, float *basis, int basisNum, int tbLen, int Nvolume)
{
    // =============== normalize coefficients  =========================
    float CMAX;
    int i,j;
    
    for (j = 0; j < basisNum; j++)
    {
        //find the maximum coefficient
        CMAX=EPSILON;
        for (i = 0; i < Nvolume; i++)
        {
            if(CMAX < coefficients[i + j*Nvolume])
            {
                CMAX = coefficients[i + j*Nvolume];
            }
        }
        
        if(CMAX>EPSILON)
        {
            // divide the coefficients by the maximum coefficient
            for (i = 0; i < Nvolume; i++)
            {
                coefficients[i + j*Nvolume] /= CMAX;
                
                if (coefficients[i + j*Nvolume]<EPSILON )
                    coefficients[i + j*Nvolume]=EPSILON;
                else if( coefficients[i + j*Nvolume]>1.0)
                    coefficients[i + j*Nvolume]=1.0;
            }
            
            // multiply the factors by the maximum coefficient
            for (i = 0; i < tbLen; i++)
            {
                basis[i + j*tbLen] *= CMAX;
                
                if (basis[i + j*tbLen]<EPSILON)
                    basis[i + j*tbLen]=EPSILON;
            }
        }
    }
    return;
}
