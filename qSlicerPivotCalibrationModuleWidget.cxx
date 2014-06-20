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
#include <QDebug>
#include <QMessageBox>
#include <QApplication>
#include <QTimer>

// SlicerQt includes
#include "qSlicerPivotCalibrationModuleWidget.h"
#include "ui_qSlicerPivotCalibrationModule.h"

#include <vtkNew.h>
#include <vtkCommand.h>

#include <vtkMRMLLinearTransformNode.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerPivotCalibrationModuleWidgetPrivate: public Ui_qSlicerPivotCalibrationModule
{
  Q_DECLARE_PUBLIC( qSlicerPivotCalibrationModuleWidget );
protected:
  qSlicerPivotCalibrationModuleWidget* const q_ptr;
public:
  qSlicerPivotCalibrationModuleWidgetPrivate( qSlicerPivotCalibrationModuleWidget& object );
  vtkSlicerPivotCalibrationLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerPivotCalibrationModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerPivotCalibrationModuleWidgetPrivate::qSlicerPivotCalibrationModuleWidgetPrivate( qSlicerPivotCalibrationModuleWidget& object )
 : q_ptr( &object )
{
}


vtkSlicerPivotCalibrationLogic* qSlicerPivotCalibrationModuleWidgetPrivate::logic() const
{
  Q_Q( const qSlicerPivotCalibrationModuleWidget );
  return vtkSlicerPivotCalibrationLogic::SafeDownCast( q->logic() );
}


//-----------------------------------------------------------------------------
// qSlicerPivotCalibrationModuleWidget methods
//-----------------------------------------------------------------------------
qSlicerPivotCalibrationModuleWidget::qSlicerPivotCalibrationModuleWidget(QWidget* _parent) : Superclass( _parent ) , d_ptr( new qSlicerPivotCalibrationModuleWidgetPrivate(*this))
{
  this->pivotDelayTimer = new QTimer();
  pivotDelayTimer->setSingleShot(false);
  pivotDelayTimer->setInterval(1000);
  this->pivotDelayCount = 5;
  
  this->pivotSamplingTimer = new QTimer();
  pivotSamplingTimer->setSingleShot(false);
  pivotSamplingTimer->setInterval(1000);
  this->pivotSamplingCount = 5;

  this->timerSetting = 5;

  this->ActivePivot = vtkSlicerPivotCalibrationLogic::NO_PIVOT;
}

//-----------------------------------------------------------------------------
qSlicerPivotCalibrationModuleWidget::~qSlicerPivotCalibrationModuleWidget()
{
  delete this->pivotDelayTimer;
  delete this->pivotSamplingTimer;
}

//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::enter()
{
  this->Superclass::enter();
}

//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::initializeObserver(vtkMRMLNode* node)
{
  Q_D(qSlicerPivotCalibrationModuleWidget);
  
  d->logic()->ObserveTransformNode( node );
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::setup()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);
  d->setupUi(this);
  

  this->Superclass::setup();
  
  connect( pivotDelayTimer, SIGNAL( timeout() ), this, SLOT( onPivotDelayTimeout() ));
  connect( pivotSamplingTimer, SIGNAL( timeout() ), this, SLOT( onPivotSamplingTimeout() ));
  
  connect( d->InputComboBox, SIGNAL( currentNodeChanged( vtkMRMLNode* ) ), this, SLOT( initializeObserver( vtkMRMLNode* ) ) );
  
  connect( d->MarkedPivotAButton, SIGNAL( clicked() ), this, SLOT( onMarkedPivotAStart() ) );
  connect( d->MarkedPivotBButton, SIGNAL( clicked() ), this, SLOT( onMarkedPivotBStart() ) );
  connect( d->UnmarkedPivotAButton, SIGNAL( clicked() ), this, SLOT( onUnmarkedPivotAStart() ) );
  connect( d->UnmarkedPivotBButton, SIGNAL( clicked() ), this, SLOT( onUnmarkedPivotBStart() ) );

  connect( d->CalibrateButton, SIGNAL( clicked() ), this, SLOT( onCalibrate() ) );
   
  connect( d->timerEdit, SIGNAL( valueChanged( double ) ), this, SLOT( setTimer( double ) ) );

  connect( d->ExtentXSpinBox, SIGNAL( valueChanged( int ) ), this, SLOT( onExtentXChanged( int ) ) );
  connect( d->ExtentYSpinBox, SIGNAL( valueChanged( int ) ), this, SLOT( onExtentYChanged( int ) ) );
  connect( d->DepthSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( onExtentYChanged( double ) ) );
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onMarkedPivotAStart()
{
  this->ActivePivot = vtkSlicerPivotCalibrationLogic::MARKED_PIVOT_A;
  this->onStartPivotPart();
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onMarkedPivotBStart()
{
  this->ActivePivot = vtkSlicerPivotCalibrationLogic::MARKED_PIVOT_B;
  this->onStartPivotPart();
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onUnmarkedPivotAStart()
{
  this->ActivePivot = vtkSlicerPivotCalibrationLogic::UNMARKED_PIVOT_A;
  this->onStartPivotPart();
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onUnmarkedPivotBStart()
{
  this->ActivePivot = vtkSlicerPivotCalibrationLogic::UNMARKED_PIVOT_B;
  this->onStartPivotPart();
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onCalibrate()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  d->logic()->GetImageToProbeTransform( d->OutputComboBox->currentNode() );
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onExtentXChanged( int newExtentX )
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  d->logic()->ExtentX = newExtentX;
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onExtentYChanged( int newExtentY )
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  d->logic()->ExtentY = newExtentY;
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onDepthChanged( double newDepth )
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  d->logic()->Depth = newDepth;
}



//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onStartPivotPart()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  this->pivotDelayCount = this->timerSetting;
  this->pivotSamplingCount = this->timerSetting;

  std::stringstream ss;
  ss << this->pivotDelayCount << " seconds until start";
  d->CountdownLabel->setText(ss.str().c_str());  
  
  pivotDelayTimer->start();
}



//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onPivotDelayTimeout()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);

  std::stringstream ss1;
  
  this->pivotDelayCount -= 1;
  ss1 << this->pivotDelayCount << " seconds until start";
  d->CountdownLabel->setText(ss1.str().c_str());

  if (this->pivotDelayCount <= 0)
  {
    std::stringstream ss2;
    this->pivotSamplingCount = this->timerSetting;
    ss2 << "Sampling time left: " << this->pivotSamplingCount;
    d->CountdownLabel->setText( ss2.str().c_str() );
    
    this->pivotDelayTimer->stop();
    d->logic()->RecordingState = true;
    this->pivotSamplingTimer->start();
  }

}

//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onPivotSamplingTimeout()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);
  
  std::stringstream ss;
  
  this->pivotSamplingCount -= 1;
  ss << "Sampling time left: " << this->pivotSamplingCount;
  d->CountdownLabel->setText( ss.str().c_str() );  
  
  if (this->pivotSamplingCount <= 0)
  {
    d->CountdownLabel->setText( "Sampling complete" );
    
    this->pivotSamplingTimer->stop();
    this->onPivotStop();
  }  
}



//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::onPivotStop()
{
  Q_D(qSlicerPivotCalibrationModuleWidget);
  
  d->logic()->RecordingState = false;
  double rmsError = d->logic()->AddPivot( this->ActivePivot );
  d->logic()->ClearSamples();

  if ( rmsError < 0 )
  {
    this->ActivePivot = vtkSlicerPivotCalibrationLogic::NO_PIVOT;
    return;
  }

  if ( this->ActivePivot == vtkSlicerPivotCalibrationLogic::MARKED_PIVOT_A )
  {
    std::stringstream ss;
    ss << "Complete! (RMS Error: " << rmsError << ")";
    d->MarkedPivotALabel->setText( ss.str().c_str() );
  }
  if ( this->ActivePivot == vtkSlicerPivotCalibrationLogic::MARKED_PIVOT_B )
  {
    std::stringstream ss;
    ss << "Complete! (RMS Error: " << rmsError << ")";
    d->MarkedPivotBLabel->setText( ss.str().c_str() );
  }
  if ( this->ActivePivot == vtkSlicerPivotCalibrationLogic::UNMARKED_PIVOT_A )
  {
    std::stringstream ss;
    ss << "Complete! (RMS Error: " << rmsError << ")";
    d->UnmarkedPivotALabel->setText( ss.str().c_str() );
  }
  if ( this->ActivePivot == vtkSlicerPivotCalibrationLogic::UNMARKED_PIVOT_B )
  {
    std::stringstream ss;
    ss << "Complete! (RMS Error: " << rmsError << ")";
    d->UnmarkedPivotBLabel->setText( ss.str().c_str() );
  }

  this->ActivePivot = vtkSlicerPivotCalibrationLogic::NO_PIVOT;
}


//-----------------------------------------------------------------------------
void qSlicerPivotCalibrationModuleWidget::setTimer(double time)
{
  this->timerSetting = (int)time;
}

