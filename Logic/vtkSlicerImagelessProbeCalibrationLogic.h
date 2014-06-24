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

// .NAME vtkSlicerImagelessProbeCalibrationLogic
// .SECTION Description
#ifndef __vtkSlicerImagelessProbeCalibrationLogic_h
#define __vtkSlicerImagelessProbeCalibrationLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkLandmarkTransform.h>

// MRML includes
#include "vtkMRMLSliceNode.h"
#include "qMRMLSliceWidget.h"

// STD includes
#include <cstdlib>

// VNL includes
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "vtkSlicerImagelessProbeCalibrationModuleLogicExport.h"


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_IMAGELESSPROBECALIBRATION_MODULE_LOGIC_EXPORT vtkSlicerImagelessProbeCalibrationLogic :
  public vtkSlicerModuleLogic
{
public:

  enum PivotEnumeration { NO_PIVOT, MARKED_SUPERIOR_PIVOT, MARKED_INFERIOR_PIVOT, UNMARKED_SUPERIOR_PIVOT, UNMARKED_INFERIOR_PIVOT, IMAGE_PLANE_SLIDING };

  static vtkSlicerImagelessProbeCalibrationLogic *New();
  vtkTypeMacro(vtkSlicerImagelessProbeCalibrationLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  void ObserveTransformNode( vtkMRMLNode* transformNode );
  void AddSample( vtkMatrix4x4* transformMatrix );

  double AddPivot( PivotEnumeration pivot );
  double ComputePivotPoint( vnl_vector< double >* pivot );
  void ComputeImageToProbeTransform();
  void ComputeImagePlaneNormalByPivots();
  void ComputeImagePlaneNormalBySliding();

  void GetImageToProbeTransform( vtkMRMLNode* outputNode );

  void ClearSamples();

  void ProcessMRMLNodesEvents( vtkObject* caller, unsigned long event, void* callData );
  
  //Move to protected, add accessors
  vnl_vector< double > MarkedSuperiorPivot;
  vnl_vector< double > MarkedInferiorPivot;
  vnl_vector< double > UnmarkedSuperiorPivot;
  vnl_vector< double > UnmarkedInferiorPivot;
  vnl_vector< double > ImagePlaneNormal;

  bool RecordingState;

  int ExtentX, ExtentY;
  double Depth;
  double RadiusOfCurvature;

  vtkSmartPointer< vtkLandmarkTransform > ImageToProbeTransform;

  
protected:
  vtkSlicerImagelessProbeCalibrationLogic();
  virtual ~vtkSlicerImagelessProbeCalibrationLogic();
  
  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  /*
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneEndImport();
  virtual void OnMRMLSceneStartClose();
  virtual void OnMRMLSceneEndClose();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
  */

  const char* ObservedTransformID;
  std::vector< vtkMatrix4x4* > CollectedTransforms;
  
private:

  vtkSlicerImagelessProbeCalibrationLogic(const vtkSlicerImagelessProbeCalibrationLogic&); // Not implemented
  void operator=(const vtkSlicerImagelessProbeCalibrationLogic&);               // Not implemented
};

#endif
