#include "lowMachFlowSolver.hpp"
#include "lowMachFlow.h"
#include "utilities/petscError.hpp"

ablate::finiteElement::LowMachFlowSolver::LowMachFlowSolver(std::string solverId, std::shared_ptr<domain::Region> region, std::shared_ptr<parameters::Parameters> options,
                                                            std::shared_ptr<parameters::Parameters> parameters, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> initialization,
                                                            std::vector<std::shared_ptr<boundaryConditions::BoundaryCondition>> boundaryConditions,
                                                            std::vector<std::shared_ptr<mathFunctions::FieldFunction>> auxiliaryFields,
                                                            std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions)
    : FiniteElementSolver(solverId, region, options, initialization, boundaryConditions, auxiliaryFields, exactSolutions), parameters(parameters) {}

void ablate::finiteElement::LowMachFlowSolver::Setup() {
    FiniteElementSolver::Setup();

    // Make sure that the fields
    if (VEL != subDomain->GetField("velocity").subId) {
        throw std::invalid_argument("The velocity field subId is expected to be " + std::to_string(VEL) + ", but found to be " + std::to_string(subDomain->GetField("velocity").subId));
    }
    if (PRES != subDomain->GetField("pressure").subId) {
        throw std::invalid_argument("The pressure field subId is expected to be " + std::to_string(PRES) + ", but found to be " + std::to_string(subDomain->GetField("pressure").subId));
    }
    if (TEMP != subDomain->GetField("temperature").subId) {
        throw std::invalid_argument("The temperature field subId is expected to be " + std::to_string(TEMP) + ", but found to be " + std::to_string(subDomain->GetField("temperature").subId));
    }

    {
        PetscObject pressure;
        MatNullSpace nullspacePres;

        DMGetField(subDomain->GetDM(), PRES, NULL, &pressure) >> checkError;
        MatNullSpaceCreate(PetscObjectComm(pressure), PETSC_TRUE, 0, NULL, &nullspacePres) >> checkError;
        PetscObjectCompose(pressure, "nullspace", (PetscObject)nullspacePres) >> checkError;
        MatNullSpaceDestroy(&nullspacePres) >> checkError;
    }

    PetscDS prob;
    DMGetDS(subDomain->GetDM(), &prob) >> checkError;

    // V, W, Q Test Function
    PetscDSSetResidual(prob, VTEST, LowMachFlow_vIntegrandTestFunction, LowMachFlow_vIntegrandTestGradientFunction) >> checkError;
    PetscDSSetResidual(prob, WTEST, LowMachFlow_wIntegrandTestFunction, LowMachFlow_wIntegrandTestGradientFunction) >> checkError;
    PetscDSSetResidual(prob, QTEST, LowMachFlow_qIntegrandTestFunction, NULL) >> checkError;

    PetscDSSetJacobian(prob, VTEST, VEL, LowMachFlow_g0_vu, LowMachFlow_g1_vu, NULL, LowMachFlow_g3_vu) >> checkError;
    PetscDSSetJacobian(prob, VTEST, PRES, NULL, NULL, LowMachFlow_g2_vp, NULL) >> checkError;
    PetscDSSetJacobian(prob, VTEST, TEMP, LowMachFlow_g0_vT, NULL, NULL, NULL) >> checkError;
    PetscDSSetJacobian(prob, QTEST, VEL, LowMachFlow_g0_qu, LowMachFlow_g1_qu, NULL, NULL) >> checkError;
    PetscDSSetJacobian(prob, QTEST, TEMP, LowMachFlow_g0_qT, LowMachFlow_g1_qT, NULL, NULL) >> checkError;
    PetscDSSetJacobian(prob, WTEST, VEL, LowMachFlow_g0_wu, NULL, NULL, NULL) >> checkError;
    PetscDSSetJacobian(prob, WTEST, TEMP, LowMachFlow_g0_wT, LowMachFlow_g1_wT, NULL, LowMachFlow_g3_wT) >> checkError;

    /* Setup constants */;
    PetscReal parameterArray[TOTAL_LOW_MACH_FLOW_PARAMETERS];
    parameters->Fill(TOTAL_LOW_MACH_FLOW_PARAMETERS, lowMachFlowParametersTypeNames, parameterArray, defaultParameters);
    PetscDSSetConstants(prob, TOTAL_LOW_MACH_FLOW_PARAMETERS, parameterArray) >> checkError;
}
static PetscErrorCode zero(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx) {
    PetscInt d;
    for (d = 0; d < Nc; ++d) u[d] = 0.0;
    return 0;
}

static PetscErrorCode constant(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx) {
    PetscInt d;
    for (d = 0; d < Nc; ++d) {
        u[d] = 1.0;
    }
    return 0;
}

static PetscErrorCode createPressureNullSpace(DM dm, PetscInt ofield, PetscInt nfield, MatNullSpace *nullSpace) {
    Vec vec;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    // determine the number of fields from PETSC
    PetscInt numFields;
    ierr = DMGetNumFields(dm, &numFields);
    CHKERRQ(ierr);

    // Project to the field specified in nfield
    std::vector<PetscErrorCode (*)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar *, void *)> funcs(numFields, zero);
    funcs[nfield] = constant;
    ierr = DMCreateGlobalVector(dm, &vec);
    CHKERRQ(ierr);

    // Get the label for this field
    DMLabel label;
    PetscObject field;
    ierr = DMGetField(dm, nfield, &label, &field);
    CHKERRQ(ierr);

    PetscInt ids[1] = {1};
    DMProjectFunctionLabel(dm, 0.0, label, 1, ids, -1, NULL, &funcs[0], NULL, INSERT_VALUES, vec);

    ierr = VecNormalize(vec, NULL);
    CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)vec, "Pressure Null Space");
    CHKERRQ(ierr);
    ierr = VecViewFromOptions(vec, NULL, "-pressure_nullspace_view");
    CHKERRQ(ierr);
    ierr = MatNullSpaceCreate(PetscObjectComm((PetscObject)dm), PETSC_FALSE, PRES, &vec, nullSpace);
    CHKERRQ(ierr);
    ierr = VecDestroy(&vec);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/* Make the discrete pressure discretely divergence free */
static PetscErrorCode removeDiscretePressureNullspaceOnTs(TS ts, ablate::finiteElement::LowMachFlowSolver &flow) {
    Vec u;
    PetscErrorCode ierr;
    DM dm;

    PetscFunctionBegin;
    ierr = TSGetDM(ts, &dm);
    CHKERRQ(ierr);
    ierr = TSGetSolution(ts, &u);
    CHKERRQ(ierr);
    try {
        flow.CompleteFlowInitialization(dm, u);
    } catch (std::exception &exp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB, exp.what());
    }

    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

void ablate::finiteElement::LowMachFlowSolver::Initialize() {
    ablate::finiteElement::FiniteElementSolver::Initialize();

    DMSetNullSpaceConstructor(subDomain->GetDM(), PRES, createPressureNullSpace) >> checkError;
    RegisterPreStep([&](TS ts, Solver &) { removeDiscretePressureNullspaceOnTs(ts, *this); });
}

void ablate::finiteElement::LowMachFlowSolver::CompleteFlowInitialization(DM dm, Vec u) {
    MatNullSpace nullsp;

    createPressureNullSpace(dm, PRES, PRES, &nullsp) >> checkError;
    MatNullSpaceRemove(nullsp, u) >> checkError;
    MatNullSpaceDestroy(&nullsp) >> checkError;
}

#include "registrar.hpp"
REGISTER(ablate::solver::Solver, ablate::finiteElement::LowMachFlowSolver, "incompressible FE flow", ARG(std::string, "id", "the name of the flow field"),
         OPT(ablate::domain::Region, "region", "the region to apply this solver.  Default is entire domain"),
         OPT(ablate::parameters::Parameters, "options", "options for the flow passed directly to PETSc"), ARG(ablate::parameters::Parameters, "parameters", "the flow field parameters"),
         ARG(std::vector<ablate::mathFunctions::FieldFunction>, "initialization", "the solution used to initialize the flow field"),
         ARG(std::vector<ablate::finiteElement::boundaryConditions::BoundaryCondition>, "boundaryConditions", "the boundary conditions for the flow field"),
         ARG(std::vector<ablate::mathFunctions::FieldFunction>, "auxFields", "enables and sets the update functions for the auxFields"),
         OPT(std::vector<ablate::mathFunctions::FieldFunction>, "exactSolution", "optional exact solutions that can be used for error calculations"));