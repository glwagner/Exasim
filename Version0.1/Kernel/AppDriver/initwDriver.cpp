#ifndef __INITWDRIVER
#define __INITWDRIVER

void InitwDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuInitw(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitw(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitw(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
