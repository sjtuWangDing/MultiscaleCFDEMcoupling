// #include "fvCFD.H"
// #include "singlePhaseTransportModel.H"

// #include "OFversion.H"
// #if defined(version30)
//   #include "turbulentTransportModel.H"
//   #include "pisoControl.H"
// #else
//   #include "turbulenceModel.H"
// #endif

// #if defined(versionv1606plus) || defined(version40)
//   #include "fvOptions.H"
// #else
//   #include "fvIOoptionList.H"
// #endif

// #include "fixedFluxPressureFvPatchScalarField.H"
// #include "../cfdemSolverMix/cfdemCloudMix.H"

// #if defined(anisotropicRotation)
//   #include "cfdemCloudRotation.H"
// #endif

// #if defined(superquadrics_flag)
//   #include "cfdemCloudRotationSuperquadric.H"
// #endif

// #include "implicitCouple.H"
// #include "clockModel.H"
// #include "smoothingModel.H"
// #include "forceModel.H"

// #include "./cfdTools/mixAdjustPhi.H"
// #include "dynamicFvMesh.H"

// int main(int argc, char *argv[]) {
//   #include "setRootCase.H"
//   #include "createTime.H"
//   #include "createDynamicFvMesh.H"

// #if defined(version30)
//   pisoControl piso(mesh);
//   #include "createTimeControls.H"
// #endif

//   // 创建所有场变量
//   #include "createFieldsIB.H"

//   #include "createFvOptions.H"
//   #include "initContinuityErrs.H"

//   #include "readGravitationalAcceleration.H"
//   // 在 cfdemSolverMix 求解器中必须使用 impCoupleModel
//   // 在 cfdemSolverMix 求解器中必须使用 mix force model
//   #include "./mixCheckGlobal.H"

//   // create cfdemCloud
// #if defined(anisotropicRotation)
//   cfdemCloudRotation particleCloud(mesh);
// #elif defined(superquadrics_flag)
//   cfdemCloudRotationSuperquadric particleCloud(mesh);
// #else
//   cfdemCloudMix particleCloud(mesh);
// #endif

//   #include "./mixCheckModelType.H"
//   word modelType = particleCloud.modelType();

//   // * * * * * * * * * * * * * * * * * Starting time loop * * * * * * * * * * * * * * * * * //
//   Info << "\nStarting time loop\n" << endl;

//   while(runTime.loop()) {

//     Info << "Time = " << runTime.timeName() << nl << endl;

//     // 设置动态加密网格
//     interFace = mag(mesh.lookupObject<volScalarField>("volumefractionNext"));
//     particleCloud.setMeshHasUpdatedFlag(mesh.update());

// #if defined(version30)
//     #include "readTimeControls.H"
//     #include "CourantNo.H"
//     #include "setDeltaT.H"
// #else
//     #include "readPISOControls.H"
//     #include "CourantNo.H"
// #endif

//     bool hasEvolved = particleCloud.evolve(voidfraction, volumefraction, Us, U, interFace);

//     if (hasEvolved) {
//       // 光滑隐式力场
//       particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
//     }

//     // 获取单位体积动量交换场
//     // 在 mixCheckImCoupleM.H 中限制了 momCoupleModel 只能为隐式交换, 所以在 getProperty 函数中只能是: "implicitCouple_index" or "mixImplicitCouple_index"
//     Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("mixImplicitCouple_index")).impMomSource();
//     Ksl.correctBoundaryConditions();

//     // #include "forceCheckIm.H"
//     {
//       Info << "cfdemSolverMix: forceCheckIm..." << endl;
//       // 单位体积动量交换系数与局部平均颗粒速度都是当前时间步的值, 流体速度则是上个时间步的值, 定义单位体积隐式力 fImp = Ksl * (Us - U)
//       volVectorField fImp(Ksl * (Us - U));
//       // 计算隐式力的值: forAll(fImp, cellI) { fImp[cellI] *= particleCloud.mesh().V()[cellI]; }
//       particleCloud.scaleWithVcell(fImp);
//       dimensionedVector fImpTotal = gSum(fImp);
//       Info << "TotalForceImp:  " << fImpTotal.value() << endl;
//       Info << "Warning, these values are based on latest Ksl and Us but prev. iteration U!" << endl;
//       Info << "cfdemSolverMix: forceCheckIm - done\n" << endl;
//     }
//     #include "solverDebugInfo.H"

//     if (particleCloud.solveFlow()) {
//       fvVectorMatrix UEqn
//       (
//           fvm::ddt(U)
//         + fvm::div(phi, U)
//         + turbulence->divDevReff(U)
//         ==
//           fvOptions(U)
//       );
//       UEqn.relax();
//       fvOptions.constrain(UEqn);

//       // momentumPredictor
// #if defined(version30)
//       if (piso.momentumPredictor())
// #else
//       if (momentumPredictor)
// #endif
//       {
//         if (modelType == "none") {
//           solve(UEqn == -fvc::grad(p));
//         } else {
//           FatalError << "cfdemSolverMixIB.C: " << __LINE__ << ": Not implement for modelType which is not none" << abort(FatalError);
//         }
//         fvOptions.correct(U);
//       }  // End of momentumPredictor

//       // PISO loop
// #if defined(version30)
//       while (piso.correct())
// #else
//       for (int corr = 0; corr < nCorr; corr++)
// #endif
//       {
//         volScalarField rUA = 1.0 / UEqn.A();
//         volVectorField HbyA = rUA * UEqn.H();
//         surfaceScalarField phiHbyA(
//           "phiHbyA",
// #ifdef version23
//           (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::interpolate(rUA) * fvc::ddtCorr(HbyA, phi)
// #else
//           (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, HbyA, phi);
// #endif
//         );
//         adjustPhi(phiHbyA, U, p);

// #if defined(version30)
//         while (piso.correctNonOrthogonal())
// #else 
//         for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
// #endif
//         {
//           fvScalarMatrix pEqn
//           (
//             fvm::laplacian(rUA, p) == fvc::div(phiHbyA) + particleCloud.ddtVoidfraction()
//           );
//           pEqn.setReference(pRefCell, pRefValue);

// #if defined(version30)
//           pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
//           if (piso.finalNonOrthogonalIter()) {
//             phi = phiHbyA - pEqn.flux();
//           }
// #else
//           if (corr == nCorr-1 && nonOrth == nNonOrthCorr) {
//   #if defined(versionExt32)
//             pEqn.solve(mesh.solutionDict().solver("pFinal"));
//   #else
//             pEqn.solve(mesh.solver("pFinal"));
//   #endif
//           } else {
//             pEqn.solve();
//           }
//           if (nonOrth == nNonOrthCorr) {
//             phi = phiHbyA - pEqn.flux();
//           }
// #endif
//         }
//         #include "continuityErrs.H"
//         U = HbyA - rUA * fvc::grad(p);
//         U.correctBoundaryConditions();
//       }  // End of piso loop
//       laminarTransport.correct();
//       turbulence->correct();
//       volScalarField volumefractionNext = mesh.lookupObject<volScalarField>("volumefractionNext");
//       Info << "cfdemSolverMix: calcVelocityCorrection..." << endl;
//       particleCloud.calcVelocityCorrection(p, U, phiIB, volumefractionNext);
//       Info << "cfdemSolverMix: calcVelocityCorrection - done" << endl;
//   #if defined(version22)
//       fvOptions.correct(U);
//   #endif
//     }  // End of particleCloud.solveFlow

//     runTime.write();
//     Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
//       << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
//   }
//   Info << "cfdemSolverMixIB.C: End\n" << endl;
//   return 0;
// }

#ifndef __mix_debug__
#define __mix_debug__ 1
#endif

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "OFversion.H"

#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif

#if __mix_debug__
    #include "./cfdemCloudMix.H"
#else
    #include "cfdemCloudIB.H"
#endif

#if defined(superquadrics_flag)
#include "cfdemCloudIBSuperquadric.H"
#endif
#include "implicitCouple.H"

#include "averagingModel.H"
#include "voidFractionModel.H"

#include "dynamicFvMesh.H"

#include "cellSet.H"

#if defined(version22)
    #include "meshToMeshNew.H"
    #include "fvIOoptionList.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif

    #include "createFieldsIB.H"

    #include "initContinuityErrs.H"

    #if defined(version22)
        #include "createFvOptions.H"
    #endif

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #if defined(superquadrics_flag)
        cfdemCloudIBSuperquadric particleCloud(mesh);
    #else
        #if __mix_debug__
            cfdemCloudMix particleCloud(mesh);
        #else
            cfdemCloudIB particleCloud(mesh);
        #endif
    #endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // //=== dyM ===================
        // interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        // particleCloud.setMeshHasUpdatedFlag(mesh.update()); //dyM

        // 设置动态加密网格
        interFace = mag(mesh.lookupObject<volScalarField>("volumefractionNext"));
        particleCloud.setMeshHasUpdatedFlag(mesh.update());

        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        Info << "- evolve()" << endl;
        // particleCloud.evolve(voidfraction, interFace);
        bool hasEvolved = particleCloud.evolve(voidfraction, volumefraction, Us, U, interFace);

        // Pressure-velocity PISO corrector
        if(particleCloud.solveFlow())
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U) //fvm::ddt(voidfraction,U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
                #if defined(version22)
                ==
                fvOptions(U)
                #endif
            );

            UEqn.relax();

            #if defined(version22)
            fvOptions.constrain(UEqn);
            #endif

            #if defined(version30)
                if (piso.momentumPredictor())
            #else
                if (momentumPredictor)
            #endif
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            #if defined(version30)
                while (piso.correct())
            #else
                for (int corr=0; corr<nCorr; corr++)
            #endif
            {
                volScalarField rUA = 1.0/UEqn.A();
                surfaceScalarField rUAf(fvc::interpolate(rUA));

                U = rUA*UEqn.H();
                #ifdef version23
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + rUAf*fvc::ddtCorr(U, phi);
                #else
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, U, phi);
                #endif
                adjustPhi(phi, U, p);

                #if defined(version22)
                fvOptions.relativeFlux(phi);
                #endif

                // Non-orthogonal pressure corrector loop
                #if defined(version30)
                    while (piso.correctNonOrthogonal())
                #else
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                #endif
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    #if defined(version30)
                        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                        if (piso.finalNonOrthogonalIter())
                            phi -= pEqn.flux();
                    #else
                        if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                            #if defined(versionExt32)
                                pEqn.solve(mesh.solutionDict().solver("pFinal"));
                            #else
                                pEqn.solve(mesh.solver("pFinal"));
                            #endif
                        else
                            pEqn.solve();

                        if (nonOrth == nNonOrthCorr)
                            phi -= pEqn.flux();
                    #endif
                }

                #include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            } 
        } //end solveFlow

        laminarTransport.correct();
        turbulence->correct();

        volScalarField volumefractionNext=mesh.lookupObject<volScalarField>("volumefractionNext");
        Info << "particleCloud.calcVelocityCorrection()..." << endl;
        particleCloud.calcVelocityCorrection(p,U,phiIB,volumefractionNext);
        Info << "particleCloud.calcVelocityCorrection() - done" << endl;

        #if defined(version22)
        fvOptions.correct(U);
        #endif

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

