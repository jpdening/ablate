#include "faceMonitor.hpp"
#include "io/interval/fixedInterval.hpp"
#include "boundarySolver/boundarySolver.hpp"
#include <iostream>
#include <sstream>

ablate::monitors::FaceMonitor::FaceMonitor(std::shared_ptr<io::interval::Interval> intervalIn): interval(intervalIn ? intervalIn : std::make_shared<io::interval::FixedInterval>()){}

PetscErrorCode ablate::monitors::FaceMonitor::MonitorFace(TS ts, PetscInt step, PetscReal crtime, Vec u, void* ctx) {
    PetscFunctionBeginUser;
    //PetscErrorCode ierr;
    auto monitor = (ablate::monitors::FaceMonitor*)ctx;

    if(monitor->interval->Check(PetscObjectComm((PetscObject)ts), step, crtime)) {


    }
    PetscFunctionReturn(0);
}

void ablate::monitors::FaceMonitor::AddField(DM& dm, const char* nameField, PetscInt numComp) {
    PetscFV fvm;
    PetscFVCreate(PetscObjectComm(PetscObject(dm)), &fvm) >> checkError;
    PetscObjectSetName((PetscObject)fvm, nameField) >> checkError;
    PetscFVSetFromOptions(fvm) >> checkError;
    PetscFVSetNumComponents(fvm, numComp) >> checkError;

    DMAddField(dm, nullptr, (PetscObject)fvm) >> checkError;
    PetscFVDestroy(&fvm);
}

void ablate::monitors::FaceMonitor::Register(std::shared_ptr<ablate::solver::Solver> solverIn) {
    //Register the solver
    ablate::monitors::Monitor::Register(solverIn);

    //Copy over the DM
    DMLabel faceLabel;
    auto bSolver = std::dynamic_pointer_cast<ablate::boundarySolver::BoundarySolver>(this->GetSolver());
    if(!bSolver) {
        throw std::invalid_argument("The face monitor can only be used with the boundary solver.");
    }
    DMGetLabel(this->GetSolver()->GetSubDomain().GetDM(), bSolver->GetFieldBoundary()->GetName().c_str(), &faceLabel);
    DMPlexFilter(this->GetSolver()->GetSubDomain().GetDM(), faceLabel, 1, &faceDM);

    //Create fields on the DM
    //Add the regression rate field
    std::string regRate = "regRate";
    AddField(faceDM, regRate.c_str(), 1);

    //Add the heat flux fields
    std::string hFlux = "hFlux";
    AddField(faceDM, hFlux.c_str(), fluxNum);

    //Create the section allowing output to HDF5
    PetscSection faceSection;
    DMGetLocalSection(faceDM, &faceSection);

    //Add the regression rate
    PetscSectionSetComponentName(faceSection, Cats::regRate, 0, regRate.c_str());

    //Add the heat flux fields
    for(PetscInt c = 0; c < fluxNum; c++) {
        std::stringstream ss;
        ss << hFlux << "_" << c;
        PetscSectionSetComponentName(faceSection, Cats::hFlux, c, ss.str().c_str());
    }

    //Name the section
    PetscObjectSetName((PetscObject)faceSection, GetId().c_str());

    //Create the global vector
    DMCreateGlobalVector(faceDM, &faceVec);
    PetscObjectSetName((PetscObject)faceVec, "faceVec");
    VecSet(faceVec, 1.0);
}

void ablate::monitors::FaceMonitor::Save(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) {

    if(sequenceNumber == 0) {
        DMView(faceDM, viewer);
    }

    DMSetOutputSequenceNumber(faceDM, sequenceNumber, time);

    VecView(faceVec, viewer);
}

void ablate::monitors::FaceMonitor::Restore(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) {

    DMSetOutputSequenceNumber(faceDM, sequenceNumber, time);
    VecLoad(faceVec, viewer);
}

#include <registrar.hpp>
REGISTER(ablate::monitors::Monitor, ablate::monitors::FaceMonitor, "Outputs data along a given boundary face",
    OPT(ablate::io::interval::Interval, "interval", "This keeps track of whether or not to output this timestep"));