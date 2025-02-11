#ifndef __FHATDRIVER
#define __FHATDRIVER

void FhatDriver(dstype *fg, dstype *xg, dstype *ug1, dstype *ug2, dstype * og1, 
     dstype * og2, dstype *wg1, dstype *wg2, dstype *uh, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension            
    Int numPoints = ngf*(f2-f1);
    Int M = numPoints * ncu;
    Int N = numPoints*ncu*nd;
    Int ntau = common.ntau;
    dstype time = common.time; 
    
    if (common.extFhat==1) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            gpuFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif                            
    }
    else {     
    // left flux
    FluxDriver(fg, xg, ug1, og1, wg1, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);
    
    // right flux
    FluxDriver(&fg[N], xg, ug2, og2, wg2, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);        
    
    // Part 1: fh = fg dot nl
    AverageFlux(fg, N, backend);    
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd, backend);            

    // Part 2: Contribution due to tau*(U-UH)
    if (common.extStab>=1) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            gpuStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif                            
    }
    else if (ntau==0) { 
        AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
    }
    else if (ntau==1) { // constant scalar  
        AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
    }
    else if (ntau==ncu) { // constant diagonal tensor
        AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints, backend);
    }
    else if (ntau==ncu*ncu) { // constant full tensor      
        AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu, backend);
    }
    else {
        printf("Stabilization option is not implemented");
        exit(-1);
    }            
    }
}

#ifdef HAVE_ENZYME
void FhatDriver(dstype *fg, dstype *dfg, dstype *xg, dstype *ug1, dstype *dug1, dstype *ug2, dstype *dug2, dstype * og1, dstype * dog1,
     dstype * og2, dstype *dog2, dstype *wg1, dstype *dwg1, dstype *wg2, dstype *dwg2, dstype *uh, dstype *duh, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension            
    Int numPoints = ngf*(f2-f1);
    Int M = numPoints * ncu;
    Int N = numPoints*ncu*nd;
    Int ntau = common.ntau;
    dstype time = common.time; 
    
        ArraySetValue(dfg, 2*numPoints*ncu, 0.0, backend);


    if (common.extFhat==1) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuFhatEnzyme(fg, dfg, xg, ug1, dug1, ug2, dug2, og1, dog1, og2, dog2, wg1, dwg1, wg2, dwg2, uh, duh, 
                nl, app.tau, app.uinf, app.physicsparam, 
                time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuFhatEnzyme(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            // gpuFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    // time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
            // gpuFhatEnzyme(fg, dfg, xg, ug1, dug1, ug2, dug2, og1, dog1, og2, dog2, wg1, dwg1, wg2, dwg2, uh, duh, 
                // nl, app.tau, app.uinf, app.physicsparam, 
                // time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
            error("Custom Fhat with AD on GPU not implemented yet");
        }
#endif                            
    }
    else {     
    // left flux
    FluxDriver(fg, dfg, xg, ug1, dug1, og1, dog1, wg1, dwg1, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);

    // right flux
    FluxDriver(&fg[N], &dfg[N], xg, ug2, dug2, og2, dog2, wg2, dwg2, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);    
                          
    // Part 1: fh = fg dot nl
    AverageFlux(fg, N, backend);    
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd, backend);    

    // dfg = dfg dot nl
    AverageFlux(dfg, N, backend);    
    AverageFluxDotNormal(dfg, nl, N, M, numPoints, nd, backend);      
    // Part 2: Contribution due to tau*(U-UH)
    if (common.extStab>=1) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuStabEnzyme(fg, dfg, xg, ug1, dug1, ug2, dug2, og1, dog1, og2, dog2, wg1, dwg1, wg2, dwg2, uh, duh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuStabEnzyme(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            // gpuStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    // time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
            // opuStabEnzyme(fg, dfg, xg, ug1, dug1, ug2, dug2, og1, dog1, og2, dog2, wg1, dwg1, wg2, dwg2, uh, duh, nl, app.tau, app.uinf, app.physicsparam, 
                    // time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
            error("Custom stabilization with AD on GPU not implemented yet");
        }
#endif                            
    }
    else if (ntau==0) { 
        AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
        AddStabilization1(dfg, dug1, dug2, app.tau, M, backend);
    }
    else if (ntau==1) { // constant scalar
        AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
        AddStabilization1(dfg, dug1, dug2, app.tau, M, backend);
    }
    else if (ntau==ncu) { // constant diagonal tensor
        AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints, backend);
        AddStabilization2(dfg, dug1, dug2, app.tau, M, numPoints, backend);
    }
    else if (ntau==ncu*ncu) { // constant full tensor 
        AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu, backend);
        AddStabilization3(dfg, dug1, dug2, app.tau, M, numPoints, ncu, backend);
    }
    else {
        printf("Stabilization option is not implemented");
        exit(-1);
    }            
    }
}
#endif

#endif