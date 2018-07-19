/* sc3dclas.c - Classical 3D Point-to-point space-Charge routine. */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#include <vector>
#include <iostream>
using namespace std ;

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include <CL/opencl.h>
#endif

#define ffac (1/(4*gpt_pi*gpt_eps0))

struct sc_info
{
    cl_platform_id platform ;           // GL platform
    cl_device_id device_id;             // compute device id 
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
    
    cl_mem input;                       // device memory used for the input array
    cl_mem output;                      // device memory used for the output array
    int prevcount ;			// Previous number of particles
} ;

static void sc_sim( double t, double *x, double *p, void *info ) ;
static void sc_ter( void *info ) ;

////////////////////////////////////////////////////////////////////////////////

// Simple compute kernel which computes the classic particle-particle interaction 
//
const char *KernelSource = "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64: enable                           \n" \
"struct cudastaticpar                                                   \n" \
"{                                                                      \n" \
"   double x ;                                                          \n" \
"   double y ;                                                          \n" \
"   double z ;                                                          \n" \
"   double q ;                                                          \n" \
"   double r2 ;                                                         \n" \
"} ;                                                                    \n" \
"                                                                       \n" \
"struct cudaforce                                                       \n" \
"{                                                                      \n" \
"   double Fx ;                                                         \n" \
"   double Fy ;                                                         \n" \
"   double Fz ;                                                         \n" \
"} ;                                                                    \n" \
"                                                                       \n" \
"__kernel void cudap2pstatic(                                           \n" \
"   __global struct cudastaticpar* input,                               \n" \
"   __global struct cudaforce* output,                                  \n" \
"   const unsigned int count,                                           \n" \
"   __local struct cudastaticpar* shared)                               \n" \
"{                                                                      \n" \

"   double posx = input[get_global_id(0)].x ;                           \n" \
"   double posy = input[get_global_id(0)].y ;                           \n" \
"   double posz = input[get_global_id(0)].z ;                           \n" \
"                                                                       \n" \
"   size_t tilesize = get_global_size(0) + get_num_groups(0) - 1 ;      \n" \
"   tilesize /= get_num_groups(0) ;                                     \n" \
"                                                                       \n" \
"   struct cudaforce force = {0,0,0} ;                                  \n" \
"                                                                       \n" \
"   for(size_t tile=0 ; tile<get_num_groups(0) ; tile++)                \n" \
"   {                                                                   \n" \
"     int offset = tile*tilesize + get_local_id(0) ;                    \n" \
"     shared[get_local_id(0)] = input[offset] ;                         \n" \
"                                                                       \n" \
"     barrier(CLK_LOCAL_MEM_FENCE) ;                                    \n" \
"                                                                       \n" \
"     for(int j=0 ; j<get_local_size(0) ; j++)                          \n" \
"     {                                                                 \n" \
"       double dx = posx-shared[j].x ;                                  \n" \
"       double dy = posy-shared[j].y ;                                  \n" \
"       double dz = posz-shared[j].z ;                                  \n" \
"       double d2 = dx*dx+dy*dy+dz*dz + shared[j].r2 ;                  \n" \
"       double d3 = d2*sqrt(d2) ;                                       \n" \
"       double fac = shared[j].q / d3 ;                                 \n" \
"       force.Fx += fac*dx ;                                            \n" \
"       force.Fy += fac*dy ;                                            \n" \
"       force.Fz += fac*dz ;                                            \n" \
"     }                                                                 \n" \
"                                                                       \n" \
"     barrier(CLK_LOCAL_MEM_FENCE) ;                                    \n" \
"   }                                                                   \n" \
"                                                                       \n" \
"   output[get_global_id(0)] = force ;                                  \n" \
"}                                                                      \n" \
"\n";

const char *KernelSourceNonTiled = "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64: enable                           \n" \
"struct cudastaticpar                                                   \n" \
"{                                                                      \n" \
"   double x ;                                                          \n" \
"   double y ;                                                          \n" \
"   double z ;                                                          \n" \
"   double q ;                                                          \n" \
"   double r2 ;                                                         \n" \
"} ;                                                                    \n" \
"                                                                       \n" \
"struct cudaforce                                                       \n" \
"{                                                                      \n" \
"   double Fx ;                                                         \n" \
"   double Fy ;                                                         \n" \
"   double Fz ;                                                         \n" \
"} ;                                                                    \n" \
"                                                                       \n" \
"__kernel void cudap2pstatic(                                           \n" \
"   __global struct cudastaticpar* input,                               \n" \
"   __global struct cudaforce* output,                                  \n" \
"   const unsigned int count)                                           \n" \
"{                                                                      \n" \
"     struct cudaforce force = {0,0,0} ;                                \n" \
"     for(int j=0 ; j<get_global_size(0) ; j++)                         \n" \
"     {                                                                 \n" \
"       double dx = input[get_global_id(0)].x-input[j].x ;              \n" \
"       double dy = input[get_global_id(0)].y-input[j].y ;              \n" \
"       double dz = input[get_global_id(0)].z-input[j].z ;              \n" \
"       double d2 = dx*dx+dy*dy+dz*dz + input[j].r2 ;                   \n" \
"       double d3 = d2*sqrt(d2) ;                                       \n" \
"       double fac = input[j].q / d3 ;                                  \n" \
"       force.Fx += fac*dx ;                                            \n" \
"       force.Fy += fac*dy ;                                            \n" \
"       force.Fz += fac*dz ;                                            \n" \
"     }                                                                 \n" \
"                                                                       \n" \
"   output[get_global_id(0)] = force ;                                  \n" \
"}                                                                      \n" \
"\n";

struct cudastaticpar
{
   cl_double x ;
   cl_double y ;
   cl_double z ;
   cl_double q ;
   cl_double r2 ;
} ;

struct cudaforce
{
   cl_double Fx ;
   cl_double Fy ;
   cl_double Fz ;
} ;

////////////////////////////////////////////////////////////////////////////////

void spacecharge3Dgpu_init(gptinit *init)
{
  struct sc_info *info ;
  int numarg ;

  numarg = gptgetargnum(init) ;
  if( numarg!=0 )
    gpterror( "Syntax: %s()", gptgetname(init) ) ;

  info = (struct sc_info *)gptmalloc(sizeof(struct sc_info)) ;

    int err;                            // error code returned from api calls
    
    // Get Platform and print info
    err = clGetPlatformIDs(1,&info->platform,NULL) ;
    if( err != CL_SUCCESS )
        gpterror("Error %d: Failed to obtain CL platform!\n", err) ;

    // Connect to a compute device and display name
    err = clGetDeviceIDs(info->platform, CL_DEVICE_TYPE_GPU, 1, &info->device_id, NULL);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to create a device group!\n", err);

    // Optional diagnostic output
    char devicename[128] ;
    clGetDeviceInfo(info->device_id, CL_DEVICE_NAME, 128, devicename, NULL) ;
    if(verbose)
        printf( "GPU device name: %s\n", devicename) ;
    
    // Create a compute context 
    info->context = clCreateContext(0, 1, &info->device_id, NULL, NULL, &err);
    if (!info->context)
        gpterror("Error %d: Failed to create a compute context!\n", err);
    
    // Create a command commands
    info->commands = clCreateCommandQueue(info->context, info->device_id, 0, &err);
    if (!info->commands)
        gpterror("Error %d: Failed to create a command commands!\n", err);
    
    // Create the compute program from the source buffer
    info->program = clCreateProgramWithSource(info->context, 1, (const char **) & KernelSource, NULL, &err);
    if (!info->program)
        gpterror("Error %d: Failed to create compute program!\n", err);
    
    // Build the program executable
    err = clBuildProgram(info->program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];
        
        fprintf(stderr,"Error %d: Failed to build program executable!\n", err);
        clGetProgramBuildInfo(info->program, info->device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        gpterror("%s\n", buffer);
    }
    
    // Create the compute kernel in the program we wish to run
    info->kernel = clCreateKernel(info->program, "cudap2pstatic", &err);
    if (!info->kernel || err != CL_SUCCESS)
        gpterror("Error %d: Failed to create compute kernel!\n", err);
    
    // Memory
    info->prevcount = 0 ;

    // Register callbacks
    odeaddfprfunction( ODEFNC_USR, sc_sim, info ) ;
    gptaddmainfunction( GPTMAINFNC_TER, sc_ter, info ) ;
}

static void sc_sim( double t, double *x, double *p, void *vinfo )
{
  struct sc_info *info = (struct sc_info *)vinfo ;

// fill array with particles
    vector<struct cudastaticpar> gpupars(numpar) ;
    for(int i=0 ; i<numpar ; i++)
    {
        gpupars[i].x = pars[i].Wr[0] ;
        gpupars[i].y = pars[i].Wr[1] ;
        gpupars[i].z = pars[i].Wr[2] ;
        gpupars[i].q = pars[i].alive ? pars[i].n * pars[i].q : 0.0 ;
        gpupars[i].r2= pars[i].r2 ;
    }

    if( numpar==0 ) return ;
    cl_int count = numpar ;

    if( info->prevcount!=numpar )
    {
        // Cleanup
        if( info->prevcount!=0 )
        {
            clReleaseMemObject(info->input);
            clReleaseMemObject(info->output);
        }

        // Create the input and output arrays in device memory for our calculation
        info->input = clCreateBuffer(info->context,  CL_MEM_READ_ONLY,  sizeof(struct cudastaticpar) * count, NULL, NULL);
        info->output = clCreateBuffer(info->context, CL_MEM_WRITE_ONLY, sizeof(struct cudaforce) * count, NULL, NULL);
        if (!info->input || !info->output)
            terminate("Error: Failed to allocate device memory!\n"); 
    }

    vector<struct cudaforce> gpuforce(numpar) ;

    // Create the input and output arrays in device memory for our calculation
    info->input = clCreateBuffer(info->context,  CL_MEM_READ_ONLY,  sizeof(struct cudastaticpar) * count, NULL, NULL);
    info->output = clCreateBuffer(info->context, CL_MEM_WRITE_ONLY, sizeof(struct cudaforce) * count, NULL, NULL);
    if (!info->input || !info->output)
        gpterror("Error: Failed to allocate device memory!\n");
    
    int err;                            // error code returned from api calls

    // Write our data set into the input array in device memory 
    err = clEnqueueWriteBuffer(info->commands, info->input, CL_TRUE, 0, sizeof(struct cudastaticpar) * count, &gpupars[0], 0, NULL, NULL);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to write to source array!\n", err);
    
    // Diagnostic output: maximum work group size for executing the kernel on the device
    size_t local;                       // local domain size for our calculation   
    err = clGetKernelWorkGroupInfo(info->kernel, info->device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to retrieve kernel work group info! %d\n");
//    if( verbose )
//       cout << "Local workgroup size:" << local << endl ;

    // Set the arguments to our compute kernel
    err = 0;
    err  = clSetKernelArg(info->kernel, 0, sizeof(cl_mem), &info->input);
    err |= clSetKernelArg(info->kernel, 1, sizeof(cl_mem), &info->output);
    err |= clSetKernelArg(info->kernel, 2, sizeof(count), &count);
    err |= clSetKernelArg(info->kernel, 3, local*sizeof(struct cudastaticpar), NULL);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to set kernel arguments!\n", err);

    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    size_t global = count ;
    err = clEnqueueNDRangeKernel(info->commands, info->kernel, 1, NULL, &global, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to execute kernel!\n", err);
    
    // Wait for the command commands to get serviced before reading back results
    err = clFinish(info->commands);
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to finish kernel!\n", err);
    
    // Read back the results from the device to verify the output
    err = clEnqueueReadBuffer( info->commands, info->output, CL_TRUE, 0, sizeof(struct cudaforce) * count, &gpuforce[0], 0, NULL, NULL );  
    if (err != CL_SUCCESS)
        gpterror("Error %d: Failed to read output array!\n", err);

    for(int i=0 ; i<numpar ; i++)
    {
        pars[i].WE[0] += ffac*gpuforce[i].Fx ;
        pars[i].WE[1] += ffac*gpuforce[i].Fy ;
        pars[i].WE[2] += ffac*gpuforce[i].Fz ;
    }
}
///////////////////////////////////////////////////////////////////////////////
// Cleanup

static void sc_ter( void *vinfo )
{
    struct sc_info *info = (struct sc_info *)vinfo ;

    clReleaseMemObject(info->input);
    clReleaseMemObject(info->output);

    clReleaseProgram(info->program);
    clReleaseKernel(info->kernel);
    clReleaseCommandQueue(info->commands);
    clReleaseContext(info->context);
}

