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

// ImagelessProbeCalibration Logic includes
#include "vtkSlicerImagelessProbeCalibrationLogic.h"

// MRML includes
#include <vtkMRMLLinearTransformNode.h>
#include "vtkMRMLScene.h"

// VTK includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkMatrix4x4.h>
#include <vtkObjectFactory.h>
#include <vtkTransform.h>
#include <vtkPoints.h>

// STD includes
#include <cassert>
#include <cmath>

// VNL includes
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/vnl_cross.h"


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerImagelessProbeCalibrationLogic);

//----------------------------------------------------------------------------
vtkSlicerImagelessProbeCalibrationLogic::vtkSlicerImagelessProbeCalibrationLogic()
{
  this->RecordingState = false;

  this->ExtentX = 512;
  this->ExtentY = 512;
  this->Depth = 50.0;

  this->MarkedSuperiorPivot = vnl_vector< double >( 3 );
  this->MarkedInferiorPivot = vnl_vector< double >( 3 );
  this->UnmarkedSuperiorPivot = vnl_vector< double >( 3 );
  this->UnmarkedInferiorPivot = vnl_vector< double >( 3 );

  this->ImageToProbeTransform = vtkSmartPointer< vtkLandmarkTransform >::New();
}

//----------------------------------------------------------------------------
vtkSlicerImagelessProbeCalibrationLogic::~vtkSlicerImagelessProbeCalibrationLogic()
{
  this->ClearSamples();
}

//----------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndImportEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  events->InsertNextValue(vtkMRMLScene::StartCloseEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::ProcessMRMLNodesEvents(vtkObject* caller, unsigned long event, void* callData)
{
  if (caller != NULL)
  {
    vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(caller);
    if ( this->RecordingState && transformNode->GetID() == this->ObservedTransformID )
    {
#ifdef TRANSFORM_NODE_MATRIX_COPY_REQUIRED
      vtkMatrix4x4* matrixCopy = vtkMatrix4x4::New();
      transformNode->GetMatrixTransformToParent(matrixCopy);
#else
      vtkMatrix4x4* matrixToParent = transformNode->GetMatrixTransformToParent();
      vtkMatrix4x4* matrixCopy = vtkMatrix4x4::New();
      matrixCopy->DeepCopy(matrixToParent);
#endif
      
      this->AddSample(matrixCopy);
    }
  }
}


//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::ObserveTransformNode(vtkMRMLNode* node)
{
  if (node != NULL)
  {
    vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(node);    
    this->ObservedTransformID = transformNode->GetID();    
    node->AddObserver(vtkMRMLLinearTransformNode::TransformModifiedEvent, (vtkCommand*) this->GetMRMLNodesCallbackCommand());
  }
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::AddSample(vtkMatrix4x4* transformMatrix)
{
  this->CollectedTransforms.push_back( transformMatrix );
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::ClearSamples()
{
  std::vector< vtkMatrix4x4* >::const_iterator it;
  for( it = this->CollectedTransforms.begin(); it != this->CollectedTransforms.end(); it++)
  {
    (*it)->Delete();
  }
  this->CollectedTransforms.clear();
}

//---------------------------------------------------------------------------
double vtkSlicerImagelessProbeCalibrationLogic::ComputePivotPoint( vnl_vector< double >* pivot )
{
  if (this->CollectedTransforms.size() == 0)
  {
    return -1;
  }

  /* Experimenting
  // Solve [ 2x 2y 2z 1 ] * [ a; b; c; d ] = [ x^2 + y^2 + z^2 ]
  // This fits the points to a sphere and calculates the centre
  vnl_matrix< double > A( this->CollectedTransforms.size(), 4 );
  vnl_vector< double > b( this->CollectedTransforms.size() );

  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    // Fill the A matrix
    A.put( i, 0, 2 * this->CollectedTransforms.at( i )->GetElement( 0, 3 ) );
    A.put( i, 1, 2 * this->CollectedTransforms.at( i )->GetElement( 1, 3 ) );
    A.put( i, 2, 2 * this->CollectedTransforms.at( i )->GetElement( 2, 3 ) );
    A.put( i, 3, 1 );

    // Fill the b vector
    b.put( i, this->CollectedTransforms.at( i )->GetElement( 0, 3 ) * this->CollectedTransforms.at( i )->GetElement( 0, 3 ) + this->CollectedTransforms.at( i )->GetElement( 1, 3 ) * this->CollectedTransforms.at( i )->GetElement( 1, 3 ) + this->CollectedTransforms.at( i )->GetElement( 2, 3 ) * this->CollectedTransforms.at( i )->GetElement( 2, 3 ) );
  }
   
    
  vnl_svd<double> pivotPointSolver( A );
  pivotPointSolver.zero_out_absolute( 1e-1 );    
  vnl_vector< double > x = pivotPointSolver.solve( b );
  */

  unsigned int rows = 3*this->CollectedTransforms.size();
  unsigned int columns = 6;

  vnl_matrix<double> A(rows, columns), minusI(3,3,0), R(3,3);
  vnl_vector<double> b(rows), x(columns), t(3);
  minusI(0, 0) = -1;
  minusI(1, 1) = -1;
  minusI(2, 2) = -1;
    
    
  std::vector<vtkMatrix4x4*>::const_iterator it, transformsEnd = this->CollectedTransforms.end();
  unsigned int currentRow;
  for(currentRow = 0, it = this->CollectedTransforms.begin(); it != transformsEnd; it++, currentRow += 3)
  {
    for (int i = 0; i < 3; i++)
    {
      t(i) = (*it)->GetElement(i, 3);
    }
    t *= -1;
    b.update(t, currentRow);

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++ )
      {
        R(i, j) = (*it)->GetElement(i, j);
      }
    }
    A.update(R, currentRow, 0);
    A.update( minusI, currentRow, 3 );    
  }
    
    
  vnl_svd<double> svdA(A);
  svdA.zero_out_absolute( 1e-1 );
  x = svdA.solve( b );

  pivot->put( 0, x.get( 0 ) );
  pivot->put( 1, x.get( 1 ) );
  pivot->put( 2, x.get( 2 ) );

  std::cout << pivot->get( 0 ) << ", " << pivot->get( 1 ) << ", " << pivot->get( 2 ) << std::endl;

  // Return RMS error
  return ( A * x - b ).rms();
}

//---------------------------------------------------------------------------
// Return true if the calibration is complete and ImageToProbe transform is available (otherwise return false)
double vtkSlicerImagelessProbeCalibrationLogic::AddPivot( PivotEnumeration pivot )
{
  if ( pivot == this->MARKED_SUPERIOR_PIVOT )
  {
    return this->ComputePivotPoint( &( this->MarkedSuperiorPivot ) );
  }
  
  if ( pivot == this->MARKED_INFERIOR_PIVOT )
  {
    return this->ComputePivotPoint( &( this->MarkedInferiorPivot ) );
  }
  
  if ( pivot == this->UNMARKED_SUPERIOR_PIVOT )
  {
    return this->ComputePivotPoint( &( this->UnmarkedSuperiorPivot ) );
  }
  
  if ( pivot == this->UNMARKED_INFERIOR_PIVOT )
  {
    return this->ComputePivotPoint( &( this->UnmarkedInferiorPivot ) );
  }





  return -1;
}

//---------------------------------------------------------------------------
// Return true if the calibration is complete and ImageToProbe transform is available (otherwise return false)
void vtkSlicerImagelessProbeCalibrationLogic::ComputeImageToProbeTransform()
{
  // First, find the image plane
  vnl_vector<double> imagePlaneNormal( 3 );
  this->ComputeImagePlaneNormal( &imagePlaneNormal );

  // Compute the marked and unmarked middle points
  vnl_vector< double > markedMiddle( 3 );
  markedMiddle.put( 0, ( this->MarkedSuperiorPivot.get( 0 ) + this->MarkedInferiorPivot.get( 0 ) ) / 2 );
  markedMiddle.put( 1, ( this->MarkedSuperiorPivot.get( 1 ) + this->MarkedInferiorPivot.get( 1 ) ) / 2 );
  markedMiddle.put( 2, ( this->MarkedSuperiorPivot.get( 2 ) + this->MarkedInferiorPivot.get( 2 ) ) / 2 );

  vnl_vector< double > unmarkedMiddle( 3 );
  unmarkedMiddle.put( 0, ( this->UnmarkedSuperiorPivot.get( 0 ) + this->UnmarkedInferiorPivot.get( 0 ) ) / 2 );
  unmarkedMiddle.put( 1, ( this->UnmarkedSuperiorPivot.get( 1 ) + this->UnmarkedInferiorPivot.get( 1 ) ) / 2 );
  unmarkedMiddle.put( 2, ( this->UnmarkedSuperiorPivot.get( 2 ) + this->UnmarkedInferiorPivot.get( 2 ) ) / 2 );

  // Compute the "far" vector by cross product
  vnl_vector< double > farVector = vnl_cross_3d( unmarkedMiddle - markedMiddle, imagePlaneNormal );
  farVector = farVector / farVector.two_norm();

  // TODO: Check the "far" vector is in the correct direction - opposite of the origin of the probe coordinate system

  vnl_vector< double > markedFar = markedMiddle + farVector * this->Depth; // DEPTH
  vnl_vector< double > unmarkedFar = unmarkedMiddle + farVector * this->Depth;

  // A properly oriented image should look like
  // X = 0, Y = 0, UN    ========== X = xmax, Y = 0, MN
  //                     ==========
  //
  //
  //
  // X = 0, Y = ymax, UF __________ X = xmax, Y = ymax, MF

  // Setup the registration
  vtkSmartPointer< vtkPoints > imagePoints = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkPoints > probePoints = vtkSmartPointer< vtkPoints >::New();

  double imagePointUN[ 3 ] = { 0, 0, 0 }; // XMAX, YMAX
  double imagePointMN[ 3 ] = { this->ExtentX, 0, 0 };
  double imagePointUF[ 3 ] = { 0, this->ExtentY, 0 };
  double imagePointMF[ 3 ] = { this->ExtentX, this->ExtentY, 0 };
  imagePoints->InsertNextPoint( imagePointUN );
  imagePoints->InsertNextPoint( imagePointMN );
  imagePoints->InsertNextPoint( imagePointUF );
  imagePoints->InsertNextPoint( imagePointMF );

  double probePointUN[ 3 ] = { unmarkedMiddle.get( 0 ), unmarkedMiddle.get( 1 ), unmarkedMiddle.get( 2 ) };
  double probePointMN[ 3 ] = { markedMiddle.get( 0 ), markedMiddle.get( 1 ), markedMiddle.get( 2 ) };
  double probePointUF[ 3 ] = { unmarkedFar.get( 0 ), unmarkedFar.get( 1 ), unmarkedFar.get( 2 ) };
  double probePointMF[ 3 ] = { markedFar.get( 0 ), markedFar.get( 1 ), markedFar.get( 2 ) };
  probePoints->InsertNextPoint( probePointUN );
  probePoints->InsertNextPoint( probePointMN );
  probePoints->InsertNextPoint( probePointUF );
  probePoints->InsertNextPoint( probePointMF );

  this->ImageToProbeTransform->SetSourceLandmarks( imagePoints );
  this->ImageToProbeTransform->SetTargetLandmarks( probePoints );
  this->ImageToProbeTransform->SetModeToSimilarity();
  this->ImageToProbeTransform->Update();

}


//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::ComputeImagePlaneNormal( vnl_vector< double >* normal )
{
  normal->put( 0, ( ( this->MarkedSuperiorPivot.get( 0 ) - this->MarkedInferiorPivot.get( 0 ) ) + ( this->UnmarkedSuperiorPivot.get( 0 ) - this->UnmarkedInferiorPivot.get( 0 ) ) ) / 2 );
  normal->put( 1, ( ( this->MarkedSuperiorPivot.get( 1 ) - this->MarkedInferiorPivot.get( 1 ) ) + ( this->UnmarkedSuperiorPivot.get( 1 ) - this->UnmarkedInferiorPivot.get( 1 ) ) ) / 2 );
  normal->put( 2, ( ( this->MarkedSuperiorPivot.get( 2 ) - this->MarkedInferiorPivot.get( 2 ) ) + ( this->UnmarkedSuperiorPivot.get( 2 ) - this->UnmarkedInferiorPivot.get( 2 ) ) ) / 2 );

  (*normal) = (*normal) / normal->two_norm();
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::GetImageToProbeTransform( vtkMRMLNode* outputNode )
{
  vtkMRMLLinearTransformNode* outputTransform = vtkMRMLLinearTransformNode::SafeDownCast( outputNode );
  if ( outputTransform == NULL )
  {
    return;
  }

  this->ComputeImageToProbeTransform();

  vtkMatrix4x4* imageToProbeMatrix = vtkMatrix4x4::New();
  this->ImageToProbeTransform->GetMatrix( imageToProbeMatrix );
  outputTransform->SetAndObserveMatrixTransformToParent( imageToProbeMatrix );
}


/*
//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}
//*/

