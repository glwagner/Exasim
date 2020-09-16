#ifndef __GETUHAT
#define __GETUHAT

int isin(Int ib, Int *a, Int n)
{
    Int in = 0;
    for (int i=0; i<n; i++)        
        if (ib == a[i])
            in = 1;
        
    return in;
}        

void UhatBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, Int ib, Int backend)
{        
    Int ncw = common.ncw;
    //Int ncq = ncu*nd;
    Int nf = f2-f1;
    Int nn = npf*nf; 
    Int nga = npf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
    Int n4 = nga*(ncx+nd+1+ncu);                // ug
    Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
    Int n6 = nga*(ncx+nd+1+ncu+nc+ncw);         // wg
    
    if (ib==0) {        
        GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0, backend);        
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);        
    }
    //else if (isin(ib, common.stgib, common.nstgib)) {        
//     else if (ib <= 2) {        
//         // synthetic turbulence generation    
// //         GetArrayAtIndex(tmp.tempg, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
// //         StgHomoTurb2(tmp.tempn, tmp.tempg, app.stgdata, app.uinf, app.stgparam, &app.stgparam[nd], 
// //                 app.physicsparam, common.time, nn, common.stgNmode, nd, backend);
//         
//         GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
//         Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
//     
//         if (nd==1) {
//             FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }
//         else if (nd==2){
//             Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
//             FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }
//         else if (nd==3) {
//             Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
//             Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
//             FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }        
//         
//         // uh
//         StgHomoTurb(&tmp.tempg[n3], &tmp.tempg[n0], app.stgdata,  app.stgparam,
//                 common.time, nn, common.stgNmode, nd, backend);
// //         StgHomoTurb2(&tmp.tempg[n3], &tmp.tempg[n0], app.stgdata, app.uinf, app.stgparam, &app.stgparam[nd], 
// //                 app.physicsparam, common.time, nn, common.stgNmode, nd, backend);
//         
//         // udg
//         GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);        
//         
//         // odg
//         if (nco>0) {
//             GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
//         }
//         if (ncw>0) {
//             GetFaceNodes(&tmp.tempg[n6], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1,backend);       
//         }
//         
//         UbouDriver(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
//                 &tmp.tempg[n6], &tmp.tempg[n1], mesh, master, app, sol, tmp, common, 
//                 npf, f1, f2, ib, backend);
//         
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
//         
// //         printArray2D(tmp.tempn,nn,ncu,backend);
// //         printArray2D(app.stgdata,common.stgNmode,10,backend);
// //         printArray2D(app.stgparam,1,2*nd,backend);
// //         printArray2D(app.uinf,1,4,backend);
// //         cout<<common.stgNmode<<endl;
// //         error("here");        
//     }
    else {
        //GetFaceNodes(tmp.tempn, sol.xdg, mesh.facecon, npf, ncx, npe, ncx, f1, f2, 1, backend);            
        //Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
        GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
    
        if (nd==1) {
            FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }
        else if (nd==3) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }        
        
        //GetElemNodes(&tmp.tempg[n3], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
        //GetFaceNodes(&tmp.tempg[n4], sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);
        GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);
        
        //if (nco>0) {
        //    GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
        //}
                
// void UbouDriver(dstype *fb, dstype *xg, dstype *udg, dstype * odg, dstype * wdg, dstype *uhg, dstype *nl, 
//         meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
//         commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)        
        UbouDriver(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n4], &tmp.tempg[n5], &tmp.tempg[n6], &tmp.tempg[n3], 
                 &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);
                       
//         if (ib==2) {
//             //printArray3D(&sol.faceg[nm+n0],ngf,nf,ncx,backend);
//             printArray3D(tmp.tempn,npf,nf,ncu,backend);
//             error("here");
//         }
        
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
    }
}

// void UhatBlock2(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
//         Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, Int ib, Int backend)
// {        
//     Int ncq = ncu*nd;
//     Int nf = f2-f1;
//     Int nn = npf*nf; 
//     Int nga = npf*nf;   
//     Int n0 = 0;                                 // xg
//     Int n1 = nga*ncx;                           // nlg
//     Int n2 = nga*(ncx+nd);                      // jac
//     Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
//     Int n4 = nga*(ncx+nd+1+ncu);                // ug
//     Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
//     Int nm = ngf*f1*(ncx+nd+1);
//     
//     if (ib==0) {        
//         GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0, backend);
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);
//     }
//     else {        
//         //GetElemNodes(&tmp.tempg[n3], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
//         
//         //GetFaceNodes(&tmp.tempg[n4], sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);
//         GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);
//         
//         //if (nco>0) {
//         //    GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
//         //}
//         
//         UbouDriver(tmp.tempn, &sol.faceg[nm+n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
//                 &sol.faceg[nm+n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);
//                         
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
//     }
// }

void GetUhat(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nbf1, Int nbf2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)   
    Int nd = common.nd;     // spatial dimension        
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face      
    
    for (Int j=nbf1; j<nbf2; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        //UhatBlock(sol, res, app, master, mesh, tmp, common, handle, f1, f2, ib, backend);
        UhatBlock(sol, res, app, master, mesh, tmp, common,
                handle, nd, npe, npf, nc, ncu, ncx, nco, f1, f2, ib, backend);
    }                           
}

#endif

