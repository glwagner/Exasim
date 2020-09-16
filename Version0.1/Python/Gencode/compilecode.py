from numpy import *
from sys import platform
import os, sys
import shutil

def compilecode(app):

    print("compile code...");

    os.chdir("app");

    codename = app['codename'];
    version = app['version'];
    appname = app['appname'];
    cpucompiler = app['cpucompiler'];
    mpicompiler = app['mpicompiler'];
    gpucompiler = app['gpucompiler'];
    cpuflags = app['cpuflags'];
    gpuflags = app['gpuflags'];

    # current directory
    cdir = os.getcwd();
    ii = cdir.find(codename)
    up = cdir[(ii+1):].count("/");
    codedir = "";
    for i in range(0,up):
        codedir = codedir + "../";

    if platform == "darwin":
        coredir = codedir + "Core/Mac/";
    elif platform == "linux" or platform == "linux2":
        coredir = codedir + "Core/Linux/";
    elif platform == "win32":
        coredir = codedir + "Core/Windows/";

    versiondir = codedir  + version;
    appdriverdir = versiondir + "/Kernel/AppDriver/";
    maindir = versiondir + "/Kernel/Main/";

    shutil.copyfile(appdriverdir + "opuApp.cpp", cdir + "/opuApp.cpp");
    shutil.copyfile(appdriverdir + "cpuApp.cpp", cdir + "/cpuApp.cpp");
    shutil.copyfile(appdriverdir + "gpuApp.cu", cdir + "/gpuApp.cu");

    compilerstr = ["" for i in range(12)]

    if  size(cpucompiler)>0:
        compilerstr[0] = cpucompiler + " -fPIC -O3 -c opuApp.cpp";
        compilerstr[1] = "ar -rvs opuApp.a opuApp.o";
    else:
        compilerstr[0] = "";
        compilerstr[1] = "";

    if  size(gpucompiler)>0:
        compilerstr[2] = gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuApp.cu";
        compilerstr[3] = "ar -rvs gpuApp.a gpuApp.o";
    else:
        compilerstr[2] = "";
        compilerstr[3] = "";

    if ( size(cpuflags)>0) and ( size(cpucompiler)>0):
        str1 = cpucompiler + " -std=c++11 " + maindir + "main.cpp " + "-o serial" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "opuCore.a " + "opuApp.a ";
        str3 = cpuflags;
        compilerstr[4] = str1 + str2 + str3;

    if ( size(cpuflags)>0) and ( size(mpicompiler)>0):
        str1 = mpicompiler + " -std=c++11 -D _MPI " + maindir + "main.cpp " + "-o mpi" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "opuCore.a " + "opuApp.a ";
        str3 = cpuflags;
        compilerstr[5] = str1 + str2 + str3;

    if ( size(cpuflags)>0) and ( size(cpucompiler)>0) and ( size(gpucompiler)>0) and ( size(gpuflags)>0):
        str1 = cpucompiler + " -std=c++11 -D _CUDA " + maindir + "main.cpp " + "-o gpu" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "gpuCore.a " + coredir + "opuCore.a opuApp.a gpuApp.a ";
        str3 = cpuflags + " " + gpuflags;
        compilerstr[6] = str1 + str2 + str3;

    if ( size(cpuflags)>0) and ( size(mpicompiler)>0) and ( size(gpucompiler)>0) and ( size(gpuflags)>0):
        str1 = mpicompiler + " -std=c++11  -D _MPI -D _CUDA " + maindir + "main.cpp " + "-o gpumpi" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "gpuCore.a " + coredir + "opuCore.a opuApp.a gpuApp.a ";
        str3 = cpuflags + " " + gpuflags;
        compilerstr[7] = str1 + str2 + str3;

    if  size(cpucompiler)>0:
        compilerstr[8] = cpucompiler + " -fPIC -O3 -c cpuApp.cpp -fopenmp";
        compilerstr[9] = "ar -rvs cpuApp.a cpuApp.o";
    else:
        compilerstr[8] = "";
        compilerstr[9] = "";

    if ( size(cpuflags)>0) and ( size(cpucompiler)>0):
        str1 = cpucompiler + " -std=c++11 " + maindir + "main.cpp" + "-o openmp" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "cpuCore.a cpuApp.a "
        str3 = "-fopenmp " + cpuflags;
        compilerstr[10] = str1 + str2 + str3;

    if ( size(cpuflags)>0) and ( size(mpicompiler)>0):
        str1 = mpicompiler + " -std=c++11 -D _MPI " + maindir + "main.cpp" + "-o openmpmpi" + appname + " ";
        str2 = coredir + "commonCore.a " + coredir + "cpuCore.a cpuApp.a ";
        str3 = "-fopenmp " + cpuflags;
        compilerstr[11] = str1 + str2 + str3;

    if app['platform'] == "cpu":
        os.system(compilerstr[0]);
        os.system(compilerstr[1]);
        if app['mpiprocs']==1:
            os.system(compilerstr[4]);
        else:
            os.system(compilerstr[5]);
    elif app['platform'] == "gpu":
        os.system(compilerstr[2]);
        os.system(compilerstr[3]);
        if app['mpiprocs']==1:
            os.system(compilerstr[6]);
        else:
            os.system(compilerstr[7]);

    os.chdir("..");

    print("Compiling done!");

    return compilerstr
