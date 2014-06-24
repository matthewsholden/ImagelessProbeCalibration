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
  this->RadiusOfCurvature = 61;
  this->Linear = false;

  this->MarkedSuperiorPivot = vnl_vector< double >( 3, 0.0 );
  this->MarkedInferiorPivot = vnl_vector< double >( 3, 0.0 );
  this->UnmarkedSuperiorPivot = vnl_vector< double >( 3, 0.0 );
  this->UnmarkedInferiorPivot = vnl_vector< double >( 3, 0.0 );
  this->ImagePlaneNormal = vnl_vector< double >( 3, 0.0 );

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
  if (this->CollectedTransforms.size() == 0)
  {
    return -1;
  }

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

  if ( pivot == this->IMAGE_PLANE_SLIDING )
  {
    this->ComputeImagePlaneNormalBySliding();
    return 0;
  }

  return -1;
}

//---------------------------------------------------------------------------
// Return true if the calibration is complete and ImageToProbe transform is available (otherwise return false)
void vtkSlicerImagelessProbeCalibrationLogic::ComputeImageToProbeTransform()
{
  // First, find the image plane
  if ( this->ImagePlaneNormal.two_norm() == 0 )
  {
    this->ComputeImagePlaneNormalByPivots();
  }

  // Compute the marked and unmarked middle points
  vnl_vector< double > markedMiddle_Probe( 3 );
  markedMiddle_Probe.put( 0, ( this->MarkedSuperiorPivot.get( 0 ) + this->MarkedInferiorPivot.get( 0 ) ) / 2 );
  markedMiddle_Probe.put( 1, ( this->MarkedSuperiorPivot.get( 1 ) + this->MarkedInferiorPivot.get( 1 ) ) / 2 );
  markedMiddle_Probe.put( 2, ( this->MarkedSuperiorPivot.get( 2 ) + this->MarkedInferiorPivot.get( 2 ) ) / 2 );

  vnl_vector< double > unmarkedMiddle_Probe( 3 );
  unmarkedMiddle_Probe.put( 0, ( this->UnmarkedSuperiorPivot.get( 0 ) + this->UnmarkedInferiorPivot.get( 0 ) ) / 2 );
  unmarkedMiddle_Probe.put( 1, ( this->UnmarkedSuperiorPivot.get( 1 ) + this->UnmarkedInferiorPivot.get( 1 ) ) / 2 );
  unmarkedMiddle_Probe.put( 2, ( this->UnmarkedSuperiorPivot.get( 2 ) + this->UnmarkedInferiorPivot.get( 2 ) ) / 2 );

  // Find the transducer face "middle"
  vnl_vector< double > faceMiddle_Probe( 3 );
  faceMiddle_Probe.put( 0, ( markedMiddle_Probe.get( 0 ) + unmarkedMiddle_Probe.get( 0 ) ) / 2 );
  faceMiddle_Probe.put( 1, ( markedMiddle_Probe.get( 1 ) + unmarkedMiddle_Probe.get( 1 ) ) / 2 );
  faceMiddle_Probe.put( 2, ( markedMiddle_Probe.get( 2 ) + unmarkedMiddle_Probe.get( 2 ) ) / 2 );

  // Compute the "far" vector by cross product
  vnl_vector< double > farVector_Probe = vnl_cross_3d( unmarkedMiddle_Probe - markedMiddle_Probe, this->ImagePlaneNormal );
  farVector_Probe = farVector_Probe / farVector_Probe.two_norm();

  vnl_vector< double > planeVector1 = farVector_Probe;
  vnl_vector< double > planeVector2 = unmarkedMiddle_Probe - markedMiddle_Probe;
  planeVector2 = planeVector2 / planeVector2.two_norm();

  // Project the middle points onto the image plane
  vnl_vector< double > markedMiddle_Plane( 2, 0.0 );
  markedMiddle_Plane.put( 0, dot_product( markedMiddle_Probe - faceMiddle_Probe, planeVector1 ) );
  markedMiddle_Plane.put( 1, dot_product( markedMiddle_Probe - faceMiddle_Probe, planeVector2 ) );

  vnl_vector< double > unmarkedMiddle_Plane( 2, 0.0 );
  unmarkedMiddle_Plane.put( 0, dot_product( unmarkedMiddle_Probe - faceMiddle_Probe, planeVector1 ) );
  unmarkedMiddle_Plane.put( 1, dot_product( unmarkedMiddle_Probe - faceMiddle_Probe, planeVector2 ) );

  vnl_vector< double > faceMiddle_Plane = ( markedMiddle_Plane + unmarkedMiddle_Plane ) / 2;

  // Find the centre of the circle
  vnl_vector< double > orthogonalVector_Plane( 2, 0.0 );
  orthogonalVector_Plane.put( 0, markedMiddle_Plane.get( 1 ) - unmarkedMiddle_Plane.get( 1 ) );
  orthogonalVector_Plane.put( 1, unmarkedMiddle_Plane.get( 0 ) - markedMiddle_Plane.get( 0 ) );
  orthogonalVector_Plane = orthogonalVector_Plane / orthogonalVector_Plane.two_norm();

  double orthogonalDistance = sqrt( this->RadiusOfCurvature * this->RadiusOfCurvature - ( markedMiddle_Plane - unmarkedMiddle_Plane ).two_norm() * ( markedMiddle_Plane - unmarkedMiddle_Plane ).two_norm() / 4 );

  vnl_vector< double > circleCentre_Plane = orthogonalVector_Plane * orthogonalDistance;

  double markedAngle_Plane = atan2( markedMiddle_Plane.get( 1 ) - circleCentre_Plane.get( 1 ), markedMiddle_Plane.get( 0 ) - circleCentre_Plane.get( 0 ) );
  double unmarkedAngle_Plane = atan2( unmarkedMiddle_Plane.get( 1 ) - circleCentre_Plane.get( 1 ), unmarkedMiddle_Plane.get( 0 ) - circleCentre_Plane.get( 0 ) );
  double totalCurveAngle = abs( markedAngle_Plane - unmarkedAngle_Plane );

  // Unproject the circle's centre
  vnl_vector< double > circleCentre_Probe = faceMiddle_Probe + circleCentre_Plane.get( 0 ) * planeVector1 + circleCentre_Plane.get( 1 ) * planeVector2;

  // Calculate the far vectors for the marked and unmarked sides
  vnl_vector< double > markedFarVector_Probe = ( markedMiddle_Probe - circleCentre_Probe );
  markedFarVector_Probe = markedFarVector_Probe / markedFarVector_Probe.two_norm();
  
  vnl_vector< double > unmarkedFarVector_Probe = ( unmarkedMiddle_Probe - circleCentre_Probe );
  unmarkedFarVector_Probe = unmarkedFarVector_Probe / unmarkedFarVector_Probe.two_norm();

  // Calculate the proportion of image that the near x and far y points lie at
  double totalXWidth_Probe = 2 * ( ( faceMiddle_Probe - markedMiddle_Probe ).two_norm() + this->Depth * cos( asin( 1.0 ) - totalCurveAngle / 2 ) );
  double totalYHeight_Probe = this->RadiusOfCurvature - ( faceMiddle_Probe - circleCentre_Probe ).two_norm() + this->Depth;

  double nearXProportion = ( unmarkedMiddle_Probe - markedMiddle_Probe ).two_norm() / ( 2 * totalXWidth_Probe );
  double farYProportion = ( this->Depth * sin( asin( 1.0 ) - totalCurveAngle / 2 ) ) / totalYHeight_Probe;

  // TODO: Check the "far" vector is in the correct direction - opposite of the origin of the probe coordinate system

  vnl_vector< double > markedFar_Probe = markedMiddle_Probe + markedFarVector_Probe * this->Depth; // DEPTH
  vnl_vector< double > unmarkedFar_Probe = unmarkedMiddle_Probe + unmarkedFarVector_Probe * this->Depth;

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

  double imagePointUN[ 3 ] = { this->ExtentX / 2 - this->ExtentX * nearXProportion, 0, 0 }; // XMAX, YMAX
  double imagePointMN[ 3 ] = { this->ExtentX / 2 + this->ExtentX * nearXProportion, 0, 0 };
  double imagePointUF[ 3 ] = { 0, this->ExtentY * farYProportion, 0 };
  double imagePointMF[ 3 ] = { this->ExtentX, this->ExtentY * farYProportion, 0 };
  imagePoints->InsertNextPoint( imagePointUN );
  imagePoints->InsertNextPoint( imagePointMN );
  imagePoints->InsertNextPoint( imagePointUF );
  imagePoints->InsertNextPoint( imagePointMF );

  double probePointUN[ 3 ] = { unmarkedMiddle_Probe.get( 0 ), unmarkedMiddle_Probe.get( 1 ), unmarkedMiddle_Probe.get( 2 ) };
  double probePointMN[ 3 ] = { markedMiddle_Probe.get( 0 ), markedMiddle_Probe.get( 1 ), markedMiddle_Probe.get( 2 ) };
  double probePointUF[ 3 ] = { unmarkedFar_Probe.get( 0 ), unmarkedFar_Probe.get( 1 ), unmarkedFar_Probe.get( 2 ) };
  double probePointMF[ 3 ] = { markedFar_Probe.get( 0 ), markedFar_Probe.get( 1 ), markedFar_Probe.get( 2 ) };
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
void vtkSlicerImagelessProbeCalibrationLogic::ComputeImagePlaneNormalByPivots()
{
  this->ImagePlaneNormal.put( 0, ( ( this->MarkedSuperiorPivot.get( 0 ) - this->MarkedInferiorPivot.get( 0 ) ) + ( this->UnmarkedSuperiorPivot.get( 0 ) - this->UnmarkedInferiorPivot.get( 0 ) ) ) / 2 );
  this->ImagePlaneNormal.put( 1, ( ( this->MarkedSuperiorPivot.get( 1 ) - this->MarkedInferiorPivot.get( 1 ) ) + ( this->UnmarkedSuperiorPivot.get( 1 ) - this->UnmarkedInferiorPivot.get( 1 ) ) ) / 2 );
  this->ImagePlaneNormal.put( 2, ( ( this->MarkedSuperiorPivot.get( 2 ) - this->MarkedInferiorPivot.get( 2 ) ) + ( this->UnmarkedSuperiorPivot.get( 2 ) - this->UnmarkedInferiorPivot.get( 2 ) ) ) / 2 );

  this->ImagePlaneNormal = this->ImagePlaneNormal / this->ImagePlaneNormal.two_norm();
}


//---------------------------------------------------------------------------
void vtkSlicerImagelessProbeCalibrationLogic::ComputeImagePlaneNormalBySliding()
{
  // Add all the points into a matrix
  vnl_matrix< double > slidingMatrix( this->CollectedTransforms.size(), 3, 0.0 );  
  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    slidingMatrix.put( i, 0, this->CollectedTransforms.at( i )->GetElement( 0, 3 ) );
    slidingMatrix.put( i, 1, this->CollectedTransforms.at( i )->GetElement( 1, 3 ) );
    slidingMatrix.put( i, 2, this->CollectedTransforms.at( i )->GetElement( 2, 3 ) );
  }

  // Calculate the mean of this
  vnl_vector< double > meanSlidingPosition( 3, 0.0 );
  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    meanSlidingPosition = meanSlidingPosition + slidingMatrix.get_row( i );
  }
  meanSlidingPosition = meanSlidingPosition / this->CollectedTransforms.size();

  // Subtract the mean from the sliding matrix
  vnl_matrix< double > zeroMeanSlidingMatrix( this->CollectedTransforms.size(), 3 );  
  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    zeroMeanSlidingMatrix.set_row( i, slidingMatrix.get_row( i ) - meanSlidingPosition );
  }

  // Calculate the covariance matrix
  vnl_matrix< double > cov( 3, 3, 0.0 );
  cov = zeroMeanSlidingMatrix.transpose() * zeroMeanSlidingMatrix;

  // Find the eigenvalues of the covariance matrix
  // The eigenvector of the smallest eigenvalue is the normal of the plane (since eigenvalue of a symmetry matrix are orthogonal)
  vnl_matrix<double> eigenvectors( 3, 3, 0.0 );
  vnl_vector<double> eigenvalues( 3, 0.0 );
  vnl_symmetric_eigensystem_compute( cov, eigenvectors, eigenvalues );

  // Note: eigenvectors are ordered in increasing eigenvalue ( 0 = smallest, end = biggest )
  vnl_vector< double > normal_Reference( 3, 0.0 );
  normal_Reference.put( 0, eigenvectors.get( 0, 0 ) );
  normal_Reference.put( 1, eigenvectors.get( 1, 0 ) );
  normal_Reference.put( 2, eigenvectors.get( 2, 0 ) );

  // Note: The last entry here is zero! This is required because it is a vector, not a point.
  double arrayNormal_Reference[ 4 ] = { normal_Reference.get( 0 ), normal_Reference.get( 1 ), normal_Reference.get( 2 ), 0 };
  double arrayNormal_Probe[ 4 ] = { 0, 0, 0, 0 };

  // Now we need to get the normal in the probe coordinate system
  // Convert the normal to probe using all ProbeToReference for the sliding and take the mean
  vnl_matrix< double > normalTest_Probe( this->CollectedTransforms.size(), 3, 0.0 );
  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    vtkSmartPointer< vtkMatrix4x4 > currentReferenceToProbeMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
    currentReferenceToProbeMatrix->DeepCopy( this->CollectedTransforms.at( i ) );
    currentReferenceToProbeMatrix->Invert();
    currentReferenceToProbeMatrix->MultiplyPoint( arrayNormal_Reference, arrayNormal_Probe );

    normalTest_Probe.put( i, 0, arrayNormal_Probe[ 0 ] );
    normalTest_Probe.put( i, 1, arrayNormal_Probe[ 1 ] );
    normalTest_Probe.put( i, 2, arrayNormal_Probe[ 2 ] );
  }

  // Calculate the mean of this
  vnl_vector< double > normal_Probe( 3, 0.0 );
  for ( int i = 0; i < this->CollectedTransforms.size(); i++ )
  {
    normal_Probe = normal_Probe + normalTest_Probe.get_row( i );
  }
  normal_Probe = normal_Probe / this->CollectedTransforms.size();

  // Set the output
  this->ImagePlaneNormal.put( 0, normal_Probe.get( 0 ) );
  this->ImagePlaneNormal.put( 1, normal_Probe.get( 1 ) );
  this->ImagePlaneNormal.put( 2, normal_Probe.get( 2 ) );
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

