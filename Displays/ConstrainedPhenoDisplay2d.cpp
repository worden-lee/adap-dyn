#include "ConstrainedPhenoDisplay2d.h"

// Constraint2dController::Constraint2dController(double period,
// 					       double filePeriod,
// 					       Site *s)
//   : DisplayController<ConstrainedPhenoDisplay2d>(period,filePeriod,s)
// {}

void Constraint2dController::createDisplay()
{ display = new ConstrainedPhenoDisplay2d(this);
  display->initialize();
}

Constraint2dController::~Constraint2dController()
{ delete display; }
