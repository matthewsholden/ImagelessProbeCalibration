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

#ifndef __qSlicerImagelessProbeCalibrationModuleWidget_h
#define __qSlicerImagelessProbeCalibrationModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"
#include "vtkSlicerImagelessProbeCalibrationLogic.h"

#include <QTimer>

#include <vtkMRMLNode.h>

#include "qSlicerImagelessProbeCalibrationModuleExport.h"

class qSlicerImagelessProbeCalibrationModuleWidgetPrivate;
class vtkAddSampleCallback;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_IMAGELESSPROBECALIBRATION_EXPORT qSlicerImagelessProbeCalibrationModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerImagelessProbeCalibrationModuleWidget(QWidget *parent=0);
  virtual ~qSlicerImagelessProbeCalibrationModuleWidget();
  
  virtual void enter();

protected slots:
  void initializeObserver(vtkMRMLNode*);
  void onStartPivotPart();
  void onPivotStop();

  void onMarkedSuperiorPivotStart();
  void onMarkedInferiorPivotStart();
  void onUnmarkedSuperiorPivotStart();
  void onUnmarkedInferiorPivotStart();

  void onCalibrate();
  
  void setTimer( double );
  void onExtentXChanged( int );
  void onExtentYChanged( int );
  void onDepthChanged( double );
  
  void onPivotDelayTimeout();
  void onPivotSamplingTimeout();
  
protected:
  QScopedPointer<qSlicerImagelessProbeCalibrationModuleWidgetPrivate> d_ptr;

  virtual void setup();
  
  int timerSetting;
  
  QTimer* pivotDelayTimer;
  int pivotDelayCount;
  
  QTimer* pivotSamplingTimer;
  int pivotSamplingCount;

  vtkSlicerImagelessProbeCalibrationLogic::PivotEnumeration ActivePivot;


private:
  Q_DECLARE_PRIVATE(qSlicerImagelessProbeCalibrationModuleWidget);
  Q_DISABLE_COPY(qSlicerImagelessProbeCalibrationModuleWidget);
};

#endif
