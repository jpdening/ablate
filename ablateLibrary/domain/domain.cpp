#include "domain.hpp"
#include <typeinfo>
#include <utilities/mpiError.hpp>
#include "monitors/logs/stdOut.hpp"
#include "solver/solver.hpp"
#include "subDomain.hpp"
#include "utilities/demangler.hpp"
#include "utilities/petscError.hpp"

ablate::domain::Domain::Domain(DM dmIn, std::string name, std::vector<std::shared_ptr<FieldDescriptor>> fieldDescriptorsIn, std::vector<std::shared_ptr<modifiers::Modifier>> modifiersIn)
    : dm(dmIn), name(name), comm(PetscObjectComm((PetscObject)dm)), fieldDescriptors(std::move(fieldDescriptorsIn)), solField(nullptr), modifiers(std::move(modifiersIn)) {
    // update the dm with the modifiers
    for (auto& modifier : modifiers) {
        modifier->Modify(dm);
    }

    // register all the solution fields with the DM, store the aux fields for later
    std::vector<std::shared_ptr<FieldDescription>> allAuxFields;
    for (const auto& fieldDescriptor : fieldDescriptors) {
        for (auto& fieldDescription : fieldDescriptor->GetFields()) {
            fieldDescription->DecompressComponents(GetDimensions());
            switch (fieldDescription->location) {
                case FieldLocation::SOL:
                    RegisterField(*fieldDescription);
                    break;
                case FieldLocation::AUX:
                    allAuxFields.push_back(fieldDescription);
                    break;
                default:
                    throw std::invalid_argument("Unknown Field Location for " + fieldDescription->name);
            }
        }
    }

    // Set up the global DS
    DMCreateDS(dm) >> checkError;

    // based upon the ds divisions in the dm, create a subDomain for each
    PetscInt numberDS;
    DMGetNumDS(dm, &numberDS) >> checkError;

    // March over each ds and create a subDomain
    for (PetscInt ds = 0; ds < numberDS; ds++) {
        subDomains.emplace_back(std::make_shared<ablate::domain::SubDomain>(*this, ds, allAuxFields));
    }
}

ablate::domain::Domain::~Domain() {
    // clean up the petsc objects
    if (solField) {
        VecDestroy(&solField) >> checkError;
    }
}

void ablate::domain::Domain::RegisterField(const ablate::domain::FieldDescription& fieldDescription) {
    // make sure that this is a solution field
    if (fieldDescription.location != FieldLocation::SOL) {
        throw std::invalid_argument("The field must be FieldLocation::SOL to be registered with the domain");
    }

    // Look up the label for this field
    DMLabel label = nullptr;
    if (fieldDescription.region) {
        DMGetLabel(dm, fieldDescription.region->GetName().c_str(), &label) >> checkError;
        if (label == nullptr) {
            throw std::invalid_argument("Cannot locate label " + fieldDescription.region->GetName() + " for field " + fieldDescription.name);
        }
    }

    // Create the field and add it with the label
    auto petscField = fieldDescription.CreatePetscField(dm);

    // add to the dm
    DMAddField(dm, label, petscField);

    // Free the petsc after being added
    PetscObjectDestroy(&petscField);

    // Record the field
    fields.push_back(Field::FromFieldDescription(fieldDescription, fields.size()));
}

PetscInt ablate::domain::Domain::GetDimensions() const {
    PetscInt dim;
    DMGetDimension(dm, &dim) >> checkError;
    return dim;
}

void ablate::domain::Domain::CreateStructures() {
    // Setup the solve with the ts
    DMPlexCreateClosureIndex(dm, NULL) >> checkError;
    DMCreateGlobalVector(dm, &(solField)) >> checkError;
    PetscObjectSetName((PetscObject)solField, "solution") >> checkError;

    // add the names to each of the components in the dm section
    PetscSection section;
    DMGetLocalSection(dm, &section) >> checkError;
    for (const auto& field : fields) {
        if (field.numberComponents > 1) {
            for (PetscInt c = 0; c < field.numberComponents; c++) {
                PetscSectionSetComponentName(section, field.id, c, field.components[c].c_str()) >> checkError;
            }
        }
    }
}

std::shared_ptr<ablate::domain::SubDomain> ablate::domain::Domain::GetSubDomain(std::shared_ptr<domain::Region> region) {
    // Check to see if there is a label for this region
    if (region) {
        // March over each ds region, and return the subdomain if this region is inside of any subDomain region
        for (const auto& subDomain : subDomains) {
            if (subDomain->InRegion(*region)) {
                return subDomain;
            }
        }
        throw std::runtime_error("Unable to locate subDomain for region " + region->ToString());

    } else {
        // Get the only subDomain
        if (subDomains.size() > 1) {
            throw std::runtime_error("More than one DS was created, the region is expected to be defined.");
        }
        return subDomains.front();
    }
}

void ablate::domain::Domain::InitializeSubDomains(std::vector<std::shared_ptr<solver::Solver>> solvers, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> initializations) {
    // determine the number of fields
    for (auto& solver : solvers) {
        solver->Register(GetSubDomain(solver->GetRegion()));
    }

    // Setup each of the fields
    for (auto& solver : solvers) {
        solver->Setup();
    }
    // Create the global structures
    CreateStructures();
    for (auto& subDomain : subDomains) {
        subDomain->CreateSubDomainStructures();
    }

    // Set the initial conditions for each field specified
    PetscInt numberFields;
    DMGetNumFields(dm, &numberFields) >> checkError;
    for (auto& initialization : initializations) {
        // Size up the field projects
        std::vector<mathFunctions::PetscFunction> fieldFunctions(numberFields, nullptr);
        std::vector<void*> fieldContexts(numberFields, nullptr);

        auto fieldId = GetField(initialization->GetName());
        fieldContexts[fieldId.id] = initialization->GetSolutionField().GetContext();
        fieldFunctions[fieldId.id] = initialization->GetSolutionField().GetPetscFunction();

        // Determine where to apply this field
        DMLabel fieldLabel = nullptr;
        PetscInt fieldValue = 0;
        if (const auto& region = initialization->GetRegion()) {
            fieldValue = region->GetValue();
            DMGetLabel(dm, region->GetName().c_str(), &fieldLabel) >> checkError;
        } else {
            PetscObject fieldTemp;
            DMGetField(dm, fieldId.id, &fieldLabel, &fieldTemp) >> checkError;
            if (fieldLabel) {
                fieldValue = 1;  // this is temporary until petsc allows fields to be defined with values beside 1
            }
        }

        // Project this field
        DMProjectFunctionLabel(dm, 0.0, fieldLabel, 1, &fieldValue, -1, nullptr, fieldFunctions.data(), fieldContexts.data(), INSERT_VALUES, solField) >> checkError;
    }

    // Initialize each solver
    for (auto& solver : solvers) {
        solver->Initialize();
    }
}