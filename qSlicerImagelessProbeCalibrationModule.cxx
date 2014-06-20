/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QtPlugin>

// ExtensionTemplate Logic includes
#include <vtkSlicerImagelessProbeCalibrationLogic.h>

// ExtensionTemplate includes
#include "qSlicerImagelessProbeCalibrationModule.h"
#include "qSlicerImagelessProbeCalibrationModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerImagelessProbeCalibrationModule, qSlicerImagelessProbeCalibrationModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerImagelessProbeCalibrationModulePrivate
{
public:
  qSlicerImagelessProbeCalibrationModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerImagelessProbeCalibrationModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerImagelessProbeCalibrationModulePrivate::qSlicerImagelessProbeCalibrationModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerImagelessProbeCalibrationModule methods

//-----------------------------------------------------------------------------
qSlicerImagelessProbeCalibrationModule::qSlicerImagelessProbeCalibrationModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerImagelessProbeCalibrationModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerImagelessProbeCalibrationModule::~qSlicerImagelessProbeCalibrationModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerImagelessProbeCalibrationModule::helpText()const
{
  return "For help on how to use this module visit: <a href='https://www.assembla.com/spaces/slicerigt'>SlicerIGT</a>";
}

//-----------------------------------------------------------------------------
QString qSlicerImagelessProbeCalibrationModule::acknowledgementText()const
{
  return "This work was was funded by Cancer Care Ontario and the Ontario Consortium for Adaptive Interventions in Radiation Oncology (OCAIRO)";
}

//-----------------------------------------------------------------------------
QStringList qSlicerImagelessProbeCalibrationModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString( "Matthew S. Holden (Queen's University), Tamas Ungi (Queen's University)" );
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerImagelessProbeCalibrationModule::icon()const
{
  return QIcon(":/Icons/ImagelessProbeCalibration.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerImagelessProbeCalibrationModule::categories() const
{
  return QStringList() << "IGT";
}

//-----------------------------------------------------------------------------
QStringList qSlicerImagelessProbeCalibrationModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerImagelessProbeCalibrationModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerImagelessProbeCalibrationModule::createWidgetRepresentation()
{
  return new qSlicerImagelessProbeCalibrationModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerImagelessProbeCalibrationModule::createLogic()
{
  return vtkSlicerImagelessProbeCalibrationLogic::New();
}
