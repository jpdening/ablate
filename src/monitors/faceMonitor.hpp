#ifndef ABLATE_LIBRARY_FACEMONITOR_HPP
#define ABLATE_LIBRARY_FACEMONITOR_HPP

#include "monitor.hpp"
#include "io/serializable.hpp"
#include <string>
#include "io/interval/interval.hpp"

namespace ablate::monitors {

class FaceMonitor : public Monitor, public io::Serializable{

    //List category offsets
    enum Cats {regRate, hFlux};

   protected:
    static void AddField(DM& dm, const char* nameField, PetscInt numComp = 1);

   private:
    const std::string name = "FaceMonitor";
    const std::shared_ptr<io::interval::Interval> interval;
    const PetscInt fluxNum = 1;
    DM faceDM;
    Vec faceVec;

    static PetscErrorCode MonitorFace(TS ts, PetscInt step, PetscReal crtime, Vec u, void* ctx);

   public:
    //Overridden monitor functions
    explicit FaceMonitor(std::shared_ptr<io::interval::Interval> intervalIn = {});
    void Register(std::shared_ptr<ablate::solver::Solver> solverIn) override;
    PetscMonitorFunction GetPetscFunction() override {return MonitorFace;}

    //Overridden serializable functions
    const std::string& GetId() const override { return name; }
    void Save(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) override;
    void Restore(PetscViewer viewer, PetscInt sequenceNumber, PetscReal time) override;


};

}
#endif //ABLATE_LIBRARY_FACEMONITOR_HPP