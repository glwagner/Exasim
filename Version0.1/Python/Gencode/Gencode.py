
import os, sys, numpy
from sympy import *
from importlib import import_module
import importlib.util
from checkcompilers import checkcompilers
from setcompilers import setcompilers
from compilecode import compilecode
from runcode import runcode
from syminit import syminit
from gencodeelem import gencodeelem
from nocodeelem import nocodeelem
from gencodeelem2 import gencodeelem2
from nocodeelem2 import nocodeelem2
from gencodeelem3 import gencodeelem3
from nocodeelem3 import nocodeelem3
from gencodeface import gencodeface
from nocodeface import nocodeface
from gencodeface2 import gencodeface2
from nocodeface2 import nocodeface2

def gencode(app):

    if os.path.isdir("app")==False:
        os.mkdir("app")
    else:
        if os.path.isfile("app/opuApp.a"):
            os.remove("app/opuApp.a")
        if os.path.isfile("app/cpuApp.a"):
            os.remove("app/cpuApp.a")
        if os.path.isfile("app/gpuApp.a"):
            os.remove("app/gpuApp.a")

    xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time = syminit(app)[0:16];

    pde = import_module(app['modelfile']);
    #spec = importlib.util.spec_from_loader(app['modelfile'], loader=None);
    #pde = importlib.util.module_from_spec(spec);
    # spec = importlib.util.spec_from_file_location('pdemodel', app['modelfile'])
    # pde = importlib.util.module_from_spec(spec);
    # spec.loader.exec_module(pde);
    #import pdemodel as pde

    nc = app['nc'];
    ncu = app['ncu'];
    u = udg[0:ncu];
    u1 = udg1[ncu:];
    u2 = udg2[ncu:];
    if nc>ncu:
        q = udg[ncu:];
        q1 = udg1[ncu:];
        q2 = udg2[ncu:];
    else:
        q = [];
        q1 = [];
        q2 = [];

    if hasattr(pde, 'flux'):
        #f = pde.flux(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem("Flux", f, xdg, udg, odg, wdg, uinf, param, time);
    else:
        sys.exit('pde.flux is not defined')
    if hasattr(pde, 'source'):
        #f = pde.source(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem("Source", f, xdg, udg, odg, wdg, uinf, param, time);
    else:
        nocodeelem("Source");
    if hasattr(pde, 'mass'):
        #f = pde.mass(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem("Tdfunc", f, xdg, udg, odg, wdg, uinf, param, time);
    else:
        nocodeelem("Tdfunc");
    if hasattr(pde, 'avfield'):
        #f = pde.avfield(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Avfield", f, xdg, udg, odg, wdg, uinf, param, time);
    else:
        nocodeelem2("Avfield");
    if hasattr(pde, 'output'):
        #f = pde.output(xdg, udg, odg, wdg, uinf, param, time);
        f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
        gencodeelem2("Output", f, xdg, udg, odg, wdg, uinf, param, time);
    else:
        nocodeelem2("Output");
    if hasattr(pde, 'fbou'):
        #f = pde.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
        f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
        f = numpy.reshape(f.flatten('F'),(app['ncu'],round(f.size/app['ncu'])));
        gencodeface("Fbou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    else:
        sys.exit("pde.fbou is not defined");
    if hasattr(pde, 'ubou'):
        #f = pde.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
        f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
        f = numpy.reshape(f.flatten('F'),(app['ncu'],round(f.size/app['ncu'])));
        gencodeface("Ubou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    else:
        sys.exit("pde.ubou is not defined");
    if hasattr(pde, 'fhat'):
        #f = pde.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface2("Fhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    else:
        nocodeface2("Fhat");
    if hasattr(pde, 'uhat'):
        #f = pde.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface2("Uhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    else:
        nocodeface2("Uhat");
    if hasattr(pde, 'stab'):
        #f = pde.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
        f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
        gencodeface2("Stab", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    else:
        nocodeface2("Stab");
    if hasattr(pde, 'initu'):
        udg = pde.initu(xdg, param, uinf);
        gencodeelem3("Initu", udg, xdg, uinf, param);
    else:
        if app['model']=="ModelW" or app['model']=="modelW":
            nocodeelem3("Initu");
        else:
            error("pde.initu is not defined");
    if hasattr(pde, 'initq'):
        udg = pde.initq(xdg, param, uinf);
        gencodeelem3("Initq", udg, xdg, uinf, param);
    else:
        nocodeelem3("Initq");
    if hasattr(pde, 'initw'):
        wdg = pde.initw(xdg, param, uinf);
        gencodeelem3("Initwdg", wdg, xdg, uinf, param);
    else:
        nocodeelem3("Initwdg");
    if hasattr(pde, 'initv'):
        odg = pde.initv(xdg, param, uinf);
        gencodeelem3("Initodg", odg, xdg, uinf, param);
    else:
        nocodeelem3("Initodg");
    if hasattr(pde, 'inituq'):
        udg = pde.inituq(xdg, param, uinf);
        gencodeelem3("Initudg", udg, xdg, uinf, param);
    else:
        if app['model']=="ModelW" or app['model']=="modelW":
            error("pde.inituq is not defined");
        else:
            nocodeelem3("Initudg");
    return 0;
